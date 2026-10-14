#include "cif.hpp"
#include "pgutils.hpp"

#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cerrno>
#include <cstdint>
#include <unordered_set>

using namespace cif;

static std::string copy_escape(const std::string& s)
{
    std::string result;
    result.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '\\': result += "\\\\"; break;
            case '\t': result += "\\t";  break;
            case '\n': result += "\\n";  break;
            case '\r': result += "\\r";  break;
            default:   result += c;      break;
        }
    }
    return result;
}

static std::string quote_literal(PGconn* conn, const std::string& s)
{
    char* escaped = PQescapeLiteral(conn, s.data(), s.size());
    if (!escaped)
        throw std::runtime_error("PQescapeLiteral failed");
    std::string result(escaped);
    PQfreemem(escaped);
    return result;
}

static std::string quote_literal(PGConnPointer& conn, const std::string& s)
{
    return quote_literal(conn.get(), s);
}

static bool missing(const std::string& s)
{
    return s.empty() || s == "." || s == "?";
}

static std::string cell_string(
    const CifBlock& b,
    std::initializer_list<const char*> tags)
{
    for (const char* tag : tags) {
        if (const auto* item = b.find_item(tag)) {
            if (!missing(item->value))
                return item->value;
        }
    }

    return {};
}
static std::string num_or_null(const std::string& raw);
static std::string int_or_null(const std::string& raw)
{
    std::string s = num_or_null(raw);

    if (s == "\\N")
        return "NULL";

    try {
        size_t n = 0;
        long long v = std::stoll(s, &n);

        if (n != s.size() ||
            v < INT32_MIN ||
            v > INT32_MAX)
        {
            return "NULL";
        }

        return std::to_string(v);
    }
    catch (...) {
        return "NULL";
    }
}

static std::string sql_number(const std::string& raw)
{
    std::string s = num_or_null(raw);
    return s == "\\N" ? "NULL" : s;
}

static std::string sql_text(
    PGConnPointer& conn,
    const std::string& s)
{
    return missing(s) ? "NULL" : quote_literal(conn, s);
}

static const CifItemLoop* find_atom_loop(
    const CifBlock& block)
{
    for (const auto& loop : block.loops) {
        if (loop.check_tag("_atom_site_label") &&
            loop.check_tag("_atom_site_fract_x") &&
            loop.check_tag("_atom_site_fract_y") &&
            loop.check_tag("_atom_site_fract_z"))
        {
            return &loop;
        }
    }

    return nullptr;
}

static std::string element_from_label(
    const std::string& label)
{
    size_t i = 0;

    while (i < label.size() &&
           std::isalpha(static_cast<unsigned char>(label[i])))
    {
        ++i;
    }

    std::string symbol = label.substr(0, i);

    if (symbol.size() > 2)
        symbol.resize(2);

    if (!symbol.empty()) {
        symbol[0] = static_cast<char>(
            std::toupper(static_cast<unsigned char>(symbol[0])));

        if (symbol.size() > 1) {
            symbol[1] = static_cast<char>(
                std::tolower(static_cast<unsigned char>(symbol[1])));
        }
    }

    return symbol;
}

static inline int cell_int(const CifBlock& block, const std::string& tag) {
    const CifItem* it = block.find_item(tag);
    return it ? it->as_int() : 0;
}

static std::vector<long long> insert_structures(
    PGConnPointer& conn,
    const std::vector<CifBlock>& batch,
    const std::string& source_file)
{
    if (batch.empty())
        return {};

    std::ostringstream sql;

    sql << "INSERT INTO structures("
           "refcode, formula, space_group, space_group_num, "
           "point_group, crystal_system, is_centrosymmetric, "
           "a, b, c, alpha, beta, gamma, "
           "volume, reduced_volume, z, "
           "measurement_temperature, "
           "r_factor_gt, wr_factor_gt, source_file"
           ") VALUES ";

    for (size_t i = 0; i < batch.size(); ++i) {
        const auto& b = batch[i];

        std::string refcode = cell_string(b, {
            "_database_code_CSD",
            "_database_code_depnum_ccdc_archive"
        });

        if (refcode.empty())
            refcode = b.name;

        std::string formula = cell_string(b, {
            "_chemical_formula_sum"
        });

        std::string sg = cell_string(b, {
            "_space_group_name_H-M_alt",
            "_symmetry_space_group_name_H-M"
        });

        std::string sg_num = int_or_null(cell_string(b, {
            "_space_group_IT_number",
            "_symmetry_Int_Tables_number"
        }));

        std::string system = cell_string(b, {
            "_space_group_crystal_system",
            "_symmetry_cell_setting"
        });

        std::string point_group = cell_string(b, {
            "_space_group_point_group_H-M",
            "_symmetry_point_group_name_H-M"
        });

        std::string volume = sql_number(
            cell_string(b, {"_cell_volume"})
        );

        std::string z = int_or_null(
            cell_string(b, {"_cell_formula_units_Z"})
        );

        std::string reduced = "NULL";

        if (i)
            sql << ',';

        sql << '('
            << sql_text(conn, refcode) << ','
            << sql_text(conn, formula) << ','
            << sql_text(conn, sg) << ','
            << sg_num << ','
            << sql_text(conn, point_group) << ','
            << sql_text(conn, system) << ','

            // Центросимметричность пока не определяется.
            << "NULL,"

            << sql_number(cell_string(
                b, {"_cell_length_a"})) << ','

            << sql_number(cell_string(
                b, {"_cell_length_b"})) << ','

            << sql_number(cell_string(
                b, {"_cell_length_c"})) << ','

            << sql_number(cell_string(
                b, {"_cell_angle_alpha"})) << ','

            << sql_number(cell_string(
                b, {"_cell_angle_beta"})) << ','

            << sql_number(cell_string(
                b, {"_cell_angle_gamma"})) << ','

            << volume << ','
            << reduced << ','
            << z << ','

            << sql_number(cell_string(b, {
                "_diffrn_ambient_temperature",
                "_cell_measurement_temperature"
            })) << ','

            << sql_number(cell_string(
                b, {"_refine_ls_R_factor_gt"})) << ','

            << sql_number(cell_string(b, {
                "_refine_ls_wR_factor_gt",
                "_refine_ls_wR_factor_ref"
            })) << ','

            << quote_literal(conn, source_file)
            << ')';
    }

    sql << " RETURNING id";

    PGresult* raw = PQexec(conn.get(), sql.str().c_str());

    if (!raw) {
        throw std::runtime_error(
            "INSERT structures: PQexec returned null"
        );
    }

    PGResultPointer res(raw, PQclear);

    if (PQresultStatus(res.get()) != PGRES_TUPLES_OK) {
        throw std::runtime_error(
            "INSERT structures failed: " +
            std::string(PQerrorMessage(conn.get()))
        );
    }

    std::vector<long long> ids;

    for (int i = 0; i < PQntuples(res.get()); ++i) {
        ids.push_back(
            std::stoll(PQgetvalue(res.get(), i, 0))
        );
    }

    if (ids.size() != batch.size()) {
        throw std::runtime_error(
            "INSERT structures: id count != row count"
        );
    }

    return ids;
}

static std::string num_or_null(const std::string& raw)
{
    if (raw.empty()) return "\\N";
    std::string s = strip_bk(raw);          // "0.1234(5)" → "0.1234"
    if (s.empty() || s == "." || s == "?") return "\\N";
    try {
        std::size_t pos = 0;
        (void)std::stod(s, &pos);
        if (pos != s.size()) return "\\N";
        return s;
    } catch (...) {
        return "\\N";
    }
}

static std::vector<std::unordered_set<std::string>> copy_atoms(
                                PGConnPointer& conn,
                                const std::vector<CifBlock>& batch,
                                const std::vector<long long>& ids)
{
    PGCopyBegin(conn, "COPY atoms("
                        "structure_id,"
                        "label,"
                        "element,"
                        "x, y, z,"
                        "occupancy,"
                        "u_iso) "
                        "FROM STDIN");
    std::vector<std::unordered_set<std::string>> inserted(batch.size());

    for (size_t k = 0; k < batch.size(); ++k) {
        const CifBlock& block = batch[k];
        const long long  id   = ids[k];

        const CifItemLoop* loop = find_atom_loop(block);
        if (!loop) continue;

        const int ix = loop->get_header_index("_atom_site_fract_x");
        const int iy = loop->get_header_index("_atom_site_fract_y");
        const int iz = loop->get_header_index("_atom_site_fract_z");
        const int il = loop->get_header_index("_atom_site_label");
        const int ie = loop->get_header_index("_atom_site_type_symbol");
        const int io = loop->get_header_index("_atom_site_occupancy");
        const int iu = loop->get_header_index("_atom_site_U_iso_or_equiv");

        for (const auto& row : loop->data) {
            auto at = [&](int idx) -> std::string {
                return (idx >= 0 && idx < (int)row.size()) ? row[idx] : "";
            };

            std::string label = at(il);
            std::string element = at(ie);

            std::string sx = at(ix);
            std::string sy = at(iy);
            std::string sz = at(iz);

            if (missing(label) ||
                missing(element) ||
                num_or_null(sx) == "\\N" ||
                num_or_null(sy) == "\\N" ||
                num_or_null(sz) == "\\N")
            {
                continue;
            }

            std::string occupancy = at(io);
            std::string u         = at(iu);

            std::ostringstream line;
            line << id << '\t'
                 << copy_escape(label)   << '\t'
                 << copy_escape(element) << '\t'
                 << num_or_null(sx) << '\t'
                 << num_or_null(sy) << '\t'
                 << num_or_null(sz) << '\t'
                 << num_or_null(occupancy) << '\t'
                 << num_or_null(u) << '\n';
            
            if (!inserted[k].insert(label).second)
                continue;
            PGCopySend(conn, line.str());
        }
    }

    PGCopyEnd(conn);
    return inserted;
}

static inline const CifItemLoop* find_aniso_loop(const CifBlock& block) {
    for (const auto& loop : block.loops) {
        if (loop.check_tag("_atom_site_aniso_U_11")) return &loop;
    }
    return nullptr;
}

static void copy_adp(PGConnPointer& conn,
                    const std::vector<CifBlock>& batch,
                    const std::vector<long long>& ids,
                    const std::vector<std::unordered_set<std::string>>& atom_labels)
{
    // Сначала проверяем: есть ли хоть в одном блоке ADP?
    bool has_any = false;
    for (const auto& b : batch)
        if (find_aniso_loop(b)) { has_any = true; break; }
    if (!has_any) return;   // нет ADP — не открываем COPY зря

    PGCopyBegin(conn,
        "COPY atom_adp(structure_id, label, u11, u22, u33, u12, u13, u23) "
        "FROM STDIN");

    for (size_t k = 0; k < batch.size(); ++k) {
        const CifBlock& block = batch[k];
        const long long  id   = ids[k];

        const CifItemLoop* loop = find_aniso_loop(block);
        if (!loop) continue;
        std::unordered_set<std::string> seen;

        const int il  = loop->get_header_index("_atom_site_aniso_label");
        const int i11 = loop->get_header_index("_atom_site_aniso_U_11");
        const int i22 = loop->get_header_index("_atom_site_aniso_U_22");
        const int i33 = loop->get_header_index("_atom_site_aniso_U_33");
        const int i12 = loop->get_header_index("_atom_site_aniso_U_12");
        const int i13 = loop->get_header_index("_atom_site_aniso_U_13");
        const int i23 = loop->get_header_index("_atom_site_aniso_U_23");

        for (const auto& row : loop->data) {
            auto at = [&](int idx) -> std::string {
                return (idx >= 0 && idx < (int)row.size()) ? row[idx] : "";
            };
            std::string label = at(il);
            if (missing(label) || !atom_labels[k].count(label) || !seen.insert(label).second)
            {
                continue;
            }

            std::ostringstream line;
            line << id << '\t'
                 << copy_escape(label)    << '\t'
                 << num_or_null(at(i11))  << '\t'
                 << num_or_null(at(i22))  << '\t'
                 << num_or_null(at(i33))  << '\t'
                 << num_or_null(at(i12))  << '\t'
                 << num_or_null(at(i13))  << '\t'
                 << num_or_null(at(i23))  << '\n';
            PGCopySend(conn, line.str());
        }
    }

    PGCopyEnd(conn);
}

static const char* USAGE_MESSAGE =
    "Usage: cif_parser <input.cif> \"<conninfo>\"\n"
    "Example: cif_parser ccdc-2023-dec.cif "
    "\"host=localhost dbname=mydb user=myuser password=mypass\"\n";

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << USAGE_MESSAGE;
        return 1;
    }

    std::string input_cif = argv[1];
    std::string conninfo  = argv[2];

    PGConnPointer conn(nullptr, PQfinish);
    PGConnect(conninfo, conn);
    if (!conn || PQstatus(conn.get()) != CONNECTION_OK) {
        std::cerr << "Connection to database failed: "
                  << (conn ? PQerrorMessage(conn.get()) : "null connection")
                  << std::endl;
        return 2;
    }

    PGExec(conn, "SET synchronous_commit = OFF");
    PGExec(conn, "SET maintenance_work_mem = '1GB'");
    PGExec(conn, "BEGIN");

    std::ifstream input_file(input_cif);
    if (!input_file) {
        std::cerr << "Failed to open input CIF file: " << input_cif << std::endl;
        return 3;
    }

    CifFile cif(input_file);

    constexpr size_t BATCH_SIZE = 1000;
    std::vector<CifBlock> batch;
    batch.reserve(BATCH_SIZE);

    long long struct_counter = 0;
    long long blocks_seen    = 0;

    auto flush_batch = [&]() {
        if (batch.empty()) return;

        std::vector<long long> ids =
        insert_structures(conn, batch, input_cif);
        auto atom_labels = copy_atoms(conn, batch, ids);
        copy_adp(conn, batch, ids, atom_labels);
        batch.clear();
        PGExec(conn, "COMMIT");
        PGExec(conn, "BEGIN");
    };

    try {
        CifBlock block;
        while (cif.read_block(block)) {
            ++blocks_seen;

            if (find_atom_loop(block)) {
                ++struct_counter;
                batch.push_back(std::move(block));
                if (batch.size() >= BATCH_SIZE)
                    flush_batch();
            }

            block.clear();

            if (blocks_seen % 1000 == 0) {
                std::cerr << "Blocks: " << blocks_seen
                          << "  Structures: " << struct_counter
                          << "\r" << std::flush;
            }
        }

        flush_batch();
        PGExec(conn, "COMMIT");

        PGExec(conn, "ANALYZE structures");
        PGExec(conn, "ANALYZE atoms");

        std::cerr << "\nDone. Structures loaded: " << struct_counter << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        try { PGExec(conn, "ROLLBACK"); } catch (...) {}
        return 4;
    }

    return 0;
}