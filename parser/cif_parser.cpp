#include "import_core.hpp"
#include "pgutils.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <unordered_set>

using namespace unncdb;
namespace {
PGResultPointer query(PGConnPointer &conn, const std::string &sql,
                      const std::vector<Value> &params = {}) {
    std::vector<const char *> values;
    for (const auto &v : params)
        values.push_back(v ? v->c_str() : nullptr);
    PGResultPointer result(PQexecParams(conn.get(), sql.c_str(), static_cast<int>(params.size()),
                                        nullptr, values.data(), nullptr, nullptr, 0),
                           PQclear);
    if (!result || (PQresultStatus(result.get()) != PGRES_TUPLES_OK &&
                    PQresultStatus(result.get()) != PGRES_COMMAND_OK))
        throw std::runtime_error("Database query failed: " +
                                 std::string(PQerrorMessage(conn.get())));
    return result;
}
std::string field(const PGResultPointer &r, int row = 0, int col = 0) {
    return PQgetvalue(r.get(), row, col);
}
std::string hex_bytes(const std::string &data) {
    static const char hex[] = "0123456789abcdef";
    std::string out = "\\x";
    out.reserve(2 + 2 * data.size());
    for (unsigned char c : data) {
        out += hex[c >> 4];
        out += hex[c & 15];
    }
    return out;
}
void append_copy_value(std::string &out, const Value &value) {
    if (!value) {
        out += "\\N";
        return;
    }
    if (value->find_first_of("\\\n\r\t") == std::string::npos) {
        out += *value;
        return;
    }
    for (char c : *value)
        switch (c) {
        case '\\':
            out += "\\\\";
            break;
        case '\n':
            out += "\\n";
            break;
        case '\r':
            out += "\\r";
            break;
        case '\t':
            out += "\\t";
            break;
        default:
            out += c;
        }
}
std::string join(const std::vector<std::string> &cols) {
    std::string s;
    for (auto &c : cols) {
        if (!s.empty())
            s += ',';
        s += c;
    }
    return s;
}
void copy_rows(PGConnPointer &conn, const std::string &table, const std::vector<std::string> &cols,
               const std::vector<Row> &rows) {
    if (rows.empty())
        return;
    PGCopyBegin(conn, "COPY " + table + "(" + join(cols) + ") FROM STDIN");
    // Buffer COPY traffic instead of calling libpq once per site.
    std::string buffer;
    buffer.reserve(1024 * 1024);
    for (const auto &row : rows) {
        for (size_t j = 0; j < cols.size(); ++j) {
            if (j)
                buffer += '\t';
            auto p = row.find(cols[j]);
            if (p == row.end())
                buffer += "\\N";
            else
                append_copy_value(buffer, p->second);
        }
        buffer += '\n';
        if (buffer.size() >= 1024 * 1024) {
            PGCopySend(conn, buffer);
            buffer.clear();
        }
    }
    if (!buffer.empty())
        PGCopySend(conn, buffer);
    PGCopyEnd(conn);
}
std::vector<std::string> reserve_ids(PGConnPointer &conn, const std::string &table,
                                     const std::string &col, size_t n) {
    std::vector<std::string> ids;
    if (!n)
        return ids;
    auto r = query(
        conn, "SELECT nextval(pg_get_serial_sequence($1,$2)) FROM generate_series(1,$3::bigint)",
        {table, col, std::to_string(n)});
    ids.reserve(n);
    for (int i = 0; i < PQntuples(r.get()); ++i)
        ids.push_back(field(r, i));
    return ids;
}
const std::vector<std::string> structure_columns = {"id",
                                                    "source_block_id",
                                                    "refcode",
                                                    "formula",
                                                    "formula_moiety",
                                                    "space_group",
                                                    "space_group_num",
                                                    "point_group",
                                                    "crystal_system",
                                                    "is_centrosymmetric",
                                                    "a",
                                                    "b",
                                                    "c",
                                                    "alpha",
                                                    "beta",
                                                    "gamma",
                                                    "volume",
                                                    "reduced_volume",
                                                    "z",
                                                    "measurement_temperature",
                                                    "r_factor_gt",
                                                    "wr_factor_gt",
                                                    "r_factor_all",
                                                    "wr_factor_ref",
                                                    "source_file",
                                                    "hall_symbol",
                                                    "space_group_num_raw",
                                                    "declared_space_group_num",
                                                    "symmetry_status",
                                                    "symmetry_method",
                                                    "symmetry_version",
                                                    "radiation_probe",
                                                    "radiation_wavelength",
                                                    "description",
                                                    "refinement_details",
                                                    "disorder_status",
                                                    "extraction_status",
                                                    "atom_rows_source",
                                                    "atom_rows_stored",
                                                    "a_su",
                                                    "b_su",
                                                    "c_su",
                                                    "alpha_su",
                                                    "beta_su",
                                                    "gamma_su"};
const std::vector<std::string> atom_columns = {"atom_id",
                                               "structure_id",
                                               "source_row",
                                               "label",
                                               "element",
                                               "x",
                                               "y",
                                               "z",
                                               "u_iso",
                                               "occupancy",
                                               "raw_type_symbol",
                                               "element_source",
                                               "x_su",
                                               "y_su",
                                               "z_su",
                                               "occupancy_raw",
                                               "occupancy_state",
                                               "adp_type",
                                               "b_iso",
                                               "disorder_assembly",
                                               "disorder_group",
                                               "refinement_flags",
                                               "calc_flag"};
struct Prepared {
    RawBlock raw;
    Parsed parsed;
};
struct Stats {
    long long seen = 0, inserted = 0, reused = 0, structures = 0, atoms = 0, errors = 0,
              warnings = 0;
};
struct Options {
    std::string path, conn;
    size_t batch_size = 1000, batch_bytes = 16 * 1024 * 1024, batch_rows = 25000, max_blocks = 0;
    bool validate = false;
};
size_t positive(const std::string &s) {
    size_t pos = 0;
    auto n = std::stoull(s, &pos);
    if (pos != s.size() || s.empty() || s[0] == '-' || n == 0)
        throw std::runtime_error("Expected positive integer: " + s);
    return n;
}
Options options(int argc, char **argv) {
    if (argc < 3)
        throw std::runtime_error(
            "Usage: cif_parser input.cif 'conninfo' [--batch-size N] [--batch-bytes N] [--batch-rows N]\n       "
            "cif_parser input.cif --validate-only [--max-blocks N]");
    Options o;
    o.path = std::filesystem::absolute(argv[1]).lexically_normal().string();
    if (std::string(argv[2]) == "--validate-only")
        o.validate = true;
    else
        o.conn = argv[2];
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if (i + 1 >= argc)
            throw std::runtime_error("Missing option value: " + arg);
        size_t value = positive(argv[++i]);
        if (arg == "--batch-size")
            o.batch_size = value;
        else if (arg == "--batch-bytes")
            o.batch_bytes = value;
        else if (arg == "--batch-rows")
            o.batch_rows = value;
        else if (arg == "--max-blocks" && o.validate)
            o.max_blocks = value;
        else
            throw std::runtime_error("Unknown option: " + arg);
    }
    return o;
}
void account(const Parsed &p, Stats &stats) {
    if (p.structure)
        ++stats.structures;
    stats.atoms += static_cast<long long>(p.atoms.size());
    for (const auto &i : p.issues) {
        if (i.severity == "error")
            ++stats.errors;
        else
            ++stats.warnings;
    }
}
void summary(const Stats &s, double seconds) {
    std::cout << "{\"blocks_seen\":" << s.seen << ",\"blocks_inserted\":" << s.inserted
              << ",\"blocks_reused\":" << s.reused << ",\"structures_extracted\":" << s.structures
              << ",\"atom_rows_extracted\":" << s.atoms << ",\"errors\":" << s.errors
              << ",\"warnings\":" << s.warnings << ",\"seconds\":" << seconds << "}\n";
}
struct StopValidation {};
} // namespace
int main(int argc, char **argv) {
    PGConnPointer conn(nullptr, PQfinish);
    std::string import_id;
    Stats stats;
    double parse_seconds = 0, write_seconds = 0;
    const auto start = std::chrono::steady_clock::now();
    try {
        auto opt = options(argc, argv);
        std::ifstream input(opt.path, std::ios::binary);
        if (!input)
            throw std::runtime_error("Cannot open input: " + opt.path);
        BlockReader reader;
        if (opt.validate) {
            try {
                reader.read(input, [&](RawBlock raw) {
                    ++stats.seen;
                    auto p = parse_block(raw);
                    account(p, stats);
                    for (const auto &issue : p.issues)
                        if (issue.severity == "error")
                            std::cerr << raw.name << ": " << issue.code << ": " << issue.message
                                      << '\n';
                    if (opt.max_blocks && static_cast<size_t>(stats.seen) >= opt.max_blocks)
                        throw StopValidation{};
                });
            } catch (const StopValidation &) {
            }
            summary(
                stats,
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
            return stats.errors ? 2 : 0;
        }
        // Content identity makes retries independent of pathname and avoids duplicate loads.
        auto source_hash = file_sha256(opt.path);
        auto bytes = std::filesystem::file_size(opt.path);
        PGConnect(opt.conn, conn);
        PGExec(conn, "SET search_path TO public");
        PGExec(conn, "SET client_encoding TO 'UTF8'");
        auto version = query(conn, "SELECT version FROM schema_migrations WHERE version=2");
        if (PQntuples(version.get()) != 1)
            throw std::runtime_error("UNNCDB schema version 2 is required");
        auto src = query(
            conn,
            "INSERT INTO cif_sources(sha256,byte_size,original_path) VALUES($1,$2::bigint,$3) ON "
            "CONFLICT(sha256) DO UPDATE SET sha256=EXCLUDED.sha256 RETURNING id,integrity_status",
            {source_hash, std::to_string(bytes), opt.path});
        std::string source_id = field(src);
        if (field(src, 0, 1) == "mismatch")
            throw std::runtime_error("Source was modified during a previous import; explicit "
                                     "verification/rebuild is required");
        auto lock =
            query(conn, "SELECT pg_try_advisory_lock(818287000000000000::bigint+$1::bigint)",
                  {source_id});
        if (field(lock) != "t")
            throw std::runtime_error("Another importer is processing this source");
        auto previous = query(conn,
                              "SELECT b.ordinal,i.importer_version FROM cif_blocks b JOIN "
                              "cif_imports i ON i.id=b.import_id WHERE b.source_id=$1::bigint",
                              {source_id});
        std::unordered_set<long long> existing;
        for (int i = 0; i < PQntuples(previous.get()); ++i) {
            if (!compatible_version(field(previous, i, 1)))
                throw std::runtime_error("Source contains extraction from a different importer "
                                         "version; explicit rebuild is required");
            existing.insert(std::stoll(field(previous, i)));
        }
        // Holding the session lock proves no earlier run of this source is active.
        query(conn,
              "UPDATE cif_imports SET status='failed',finished_at=now(),error_message='Interrupted "
              "before completion; resumed by a later run' WHERE source_id=$1::bigint AND "
              "status='running'",
              {source_id});
        auto run = query(
            conn,
            "INSERT INTO cif_imports(source_id,input_path,importer_version,parser_version,status) "
            "VALUES($1::bigint,$2,$3,'gemmi-0.7.5','running') RETURNING id",
            {source_id, opt.path, VERSION});
        import_id = field(run);
        std::vector<Prepared> batch;
        size_t batch_bytes = 0, batch_rows = 0;
        auto flush = [&]() {
            if (batch.empty())
                return;
            const auto write_start = std::chrono::steady_clock::now();
            PGExec(conn, "BEGIN");
            try {
                size_t structures_count = 0, atom_count = 0;
                for (const auto &p : batch)
                    if (p.parsed.structure) {
                        ++structures_count;
                        atom_count += p.parsed.atoms.size();
                    }
                auto block_ids = reserve_ids(conn, "cif_blocks", "id", batch.size());
                auto structure_ids = reserve_ids(conn, "structures", "id", structures_count);
                auto atom_ids = reserve_ids(conn, "atoms", "atom_id", atom_count);
                std::vector<Row> blocks, structures, atoms, adp, ops, issues;
                blocks.reserve(batch.size());
                structures.reserve(structures_count);
                atoms.reserve(atom_count);
                size_t si = 0, ai = 0;
                for (size_t i = 0; i < batch.size(); ++i) {
                    auto &prepared = batch[i];
                    auto &raw = prepared.raw;
                    auto &p = prepared.parsed;
                    auto bid = block_ids[i];
                    blocks.push_back(
                        {{"id", bid},
                         {"source_id", source_id},
                         {"import_id", import_id},
                         {"ordinal", std::to_string(raw.ordinal)},
                         {"block_name", valid_utf8_text(raw.name)
                                            ? raw.name
                                            : "invalid-utf8-name-" + hex_bytes(raw.name).substr(2)},
                         {"byte_offset", std::to_string(raw.offset)},
                         {"byte_length", std::to_string(raw.raw.size())},
                         {"sha256", sha256(raw.raw)},
                         {"raw_cif", hex_bytes(raw.raw)},
                         {"parse_status", p.parse_status}});
                    for (const auto &issue : p.issues)
                        issues.push_back({{"block_id", bid},
                                          {"severity", issue.severity},
                                          {"code", issue.code},
                                          {"source_row", issue.row},
                                          {"tag", issue.tag},
                                          {"message", issue.message}});
                    if (!p.structure)
                        continue;
                    auto sid = structure_ids[si++];
                    Row s = std::move(p.fields);
                    s["id"] = sid;
                    s["source_block_id"] = bid;
                    s["source_file"] = opt.path;
                    // An extraction exception may leave only a partial set of derived fields.
                    if (!s["disorder_status"])
                        s["disorder_status"] = "unknown";
                    if (!s["extraction_status"])
                        s["extraction_status"] = "issues";
                    structures.push_back(std::move(s));
                    std::map<std::string, std::vector<std::string>> label_ids;
                    for (auto &a : p.atoms) {
                        auto id = atom_ids[ai++];
                        a["atom_id"] = id;
                        a["structure_id"] = sid;
                        if (a["label"])
                            label_ids[*a["label"]].push_back(id);
                        atoms.push_back(std::move(a));
                    }
                    for (auto &a : p.adp) {
                        a["structure_id"] = sid;
                        if (a["label"] && label_ids[*a["label"]].size() == 1)
                            a["atom_id"] = label_ids[*a["label"]][0];
                        adp.push_back(std::move(a));
                    }
                    for (auto &op : p.operations) {
                        op["structure_id"] = sid;
                        ops.push_back(std::move(op));
                    }
                }
                copy_rows(conn, "cif_blocks",
                          {"id", "source_id", "import_id", "ordinal", "block_name", "byte_offset",
                           "byte_length", "sha256", "raw_cif", "parse_status"},
                          blocks);
                copy_rows(conn, "structures", structure_columns, structures);
                copy_rows(conn, "atoms", atom_columns, atoms);
                copy_rows(conn, "atom_adp",
                          {"structure_id", "atom_id", "source_row", "label", "representation",
                           "u11", "u22", "u33", "u12", "u13", "u23", "b11", "b22", "b33", "b12",
                           "b13", "b23"},
                          adp);
                copy_rows(conn, "symmetry_operations",
                          {"structure_id", "source_row", "operation_xyz"}, ops);
                copy_rows(conn, "import_issues",
                          {"block_id", "severity", "code", "source_row", "tag", "message"}, issues);
                query(conn, "UPDATE cif_sources SET preamble=$2::bytea WHERE id=$1::bigint",
                      {source_id, hex_bytes(reader.preamble)});
                query(conn,
                      "UPDATE cif_imports SET "
                      "blocks_seen=$2::bigint,blocks_inserted=$3::bigint,blocks_reused=$4::bigint "
                      "WHERE id=$1::bigint",
                      {import_id, std::to_string(stats.seen),
                       std::to_string(stats.inserted + static_cast<long long>(batch.size())),
                       std::to_string(stats.reused)});
                PGExec(conn, "COMMIT");
                stats.inserted += static_cast<long long>(batch.size());
                batch.clear();
                batch_bytes = 0;
                batch_rows = 0;
                write_seconds +=
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - write_start)
                        .count();
                std::cerr << "Blocks: " << stats.seen << "; stored: " << stats.inserted
                          << "; reused: " << stats.reused << "; parse_s: " << std::fixed
                          << std::setprecision(1) << parse_seconds << "; write_s: " << write_seconds
                          << "; elapsed_s: "
                          << std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
                                 .count()
                          << '\n';
            } catch (...) {
                try {
                    PGExec(conn, "ROLLBACK");
                } catch (...) {
                }
                throw;
            }
        };
        Sha256 read_hash;
        bool hashed_preamble = false;
        reader.read(input, [&](RawBlock raw) {
            if (!hashed_preamble) {
                read_hash.update(reader.preamble);
                hashed_preamble = true;
            }
            read_hash.update(raw.raw);
            ++stats.seen;
            if (existing.count(raw.ordinal)) {
                ++stats.reused;
                return;
            }
            auto parse_start = std::chrono::steady_clock::now();
            auto parsed = parse_block(raw);
            parse_seconds +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - parse_start)
                    .count();
            account(parsed, stats);
            batch_rows += parsed.atoms.size() + parsed.adp.size() + parsed.operations.size() +
                          parsed.issues.size();
            batch_bytes += raw.raw.size();
            batch.push_back({std::move(raw), std::move(parsed)});
            if (batch.size() >= opt.batch_size || batch_bytes >= opt.batch_bytes ||
                batch_rows >= opt.batch_rows)
                flush();
        });
        if (!hashed_preamble)
            read_hash.update(reader.preamble);
        if (read_hash.finish() != source_hash) {
            query(conn, "UPDATE cif_sources SET integrity_status='mismatch' WHERE id=$1::bigint",
                  {source_id});
            throw std::runtime_error(
                "Source file changed during import; source integrity is marked mismatch");
        }
        flush();
        auto has_issues = query(conn,
                                "SELECT EXISTS(SELECT 1 FROM import_issues i JOIN cif_blocks b ON "
                                "b.id=i.block_id WHERE b.source_id=$1::bigint)",
                                {source_id});
        PGExec(conn, "BEGIN");
        query(conn, "UPDATE cif_sources SET preamble=$2::bytea WHERE id=$1::bigint",
              {source_id, hex_bytes(reader.preamble)});
        query(conn,
              "UPDATE cif_imports SET "
              "status=$2,finished_at=now(),blocks_seen=$3::bigint,blocks_inserted=$4::bigint,"
              "blocks_reused=$5::bigint WHERE id=$1::bigint",
              {import_id, field(has_issues) == "t" ? "complete_with_issues" : "complete",
               std::to_string(stats.seen), std::to_string(stats.inserted),
               std::to_string(stats.reused)});
        query(conn, "UPDATE cif_sources SET integrity_status='verified' WHERE id=$1::bigint",
              {source_id});
        PGExec(conn, "COMMIT");
        summary(stats,
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
        return 0;
    } catch (const std::exception &e) {
        if (conn && !import_id.empty())
            try {
                PGExec(conn, "ROLLBACK");
                query(conn,
                      "UPDATE cif_imports SET status='failed',finished_at=now(),error_message=$2 "
                      "WHERE id=$1::bigint",
                      {import_id, e.what()});
            } catch (...) {
            }
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
