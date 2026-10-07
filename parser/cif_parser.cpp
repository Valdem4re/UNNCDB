#include "cif.hpp"
#include "pgutils.hpp"

#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

static inline double cell_double(const CifBlock& block, const std::string& tag) {
    const CifItem* it = block.find_item(tag);
    return it ? it->as_double() : 0.0;
}

static inline std::string cell_string(const CifBlock& block, const std::string& tag) {
    const CifItem* it = block.find_item(tag);
    return it ? it->value : "";
}

static inline const CifItemLoop* find_atom_loop(const CifBlock& block) {
    for (const auto& loop : block.loops) {
        if (loop.check_tag("_atom_site_fract_x")) {
            return &loop;
        }
    }
    return nullptr;
}

static inline std::string element_from_label(const std::string& label) {
    size_t k = 0;
    while (k < label.size() && std::isalpha((unsigned char)label[k])) ++k;
    return label.substr(0, k);
}

static std::vector<long long> insert_structures(
        PGConnPointer& conn,
        const std::vector<CifBlock>& batch,
        const std::string& source_file)
{
    if (batch.empty()) return {};

    std::ostringstream sql;
    sql << "INSERT INTO structures("
           "refcode, formula, space_group, "
           "a, b, c, alpha, beta, gamma, volume, z, source_file"
           ") VALUES ";

    for (size_t i = 0; i < batch.size(); ++i) {
        const auto& b = batch[i];

        std::string refcode = cell_string(b, "_database_code_CSD");
        if (refcode.empty())
            refcode = cell_string(b, "_database_code_depnum_ccdc_archive");
        if (refcode.empty())
            refcode = b.name;

        std::string formula     = cell_string(b, "_chemical_formula_sum");
        std::string space_group = cell_string(b, "_symmetry_space_group_name_H-M");

        double a     = cell_double(b, "_cell_length_a");
        double bb    = cell_double(b, "_cell_length_b");
        double c     = cell_double(b, "_cell_length_c");
        double alpha = cell_double(b, "_cell_angle_alpha");
        double beta  = cell_double(b, "_cell_angle_beta");
        double gamma = cell_double(b, "_cell_angle_gamma");
        double vol   = cell_double(b, "_cell_volume");

        int z = 0;
        if (const CifItem* it = b.find_item("_cell_formula_units_Z"))
            z = it->as_int();

        if (i) sql << ',';
        sql << '('
            << quote_literal(conn, refcode)     << ','
            << quote_literal(conn, formula)     << ','
            << quote_literal(conn, space_group) << ','
            << a     << ',' << bb    << ',' << c     << ','
            << alpha << ',' << beta  << ',' << gamma << ','
            << vol   << ','
            << z     << ','
            << quote_literal(conn, source_file)
            << ')';
    }

    sql << " RETURNING id";

    PGresult* res = PQexec(conn.get(), sql.str().c_str());
    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        std::string err = PQerrorMessage(conn.get());
        PQclear(res);
        throw std::runtime_error("INSERT structures failed: " + err);
    }

    std::vector<long long> ids;
    ids.reserve(batch.size());
    const int n = PQntuples(res);
    for (int i = 0; i < n; ++i)
        ids.push_back(std::stoll(PQgetvalue(res, i, 0)));
    PQclear(res);

    if (ids.size() != batch.size())
        throw std::runtime_error(
            "INSERT structures: id count != row count");

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

static void copy_atoms(PGConnPointer& conn,
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

            std::string label   = at(il);
            std::string element = at(ie);
            if (element.empty() && !label.empty())
                element = element_from_label(label);
            if (element.empty())
                element = "X";   // неизвестный элемент — чтобы NOT NULL не упал

            std::string sx = at(ix), sy = at(iy), sz = at(iz);
            if (sx.empty() || sy.empty() || sz.empty()) continue;

            std::string occupancy = at(io);
            std::string u         = at(iu);

            std::ostringstream line;
            line << id << '\t'
                 << copy_escape(label)   << '\t'
                 << copy_escape(element) << '\t'
                 << strip_bk(sx) << '\t'
                 << strip_bk(sy) << '\t'
                 << strip_bk(sz) << '\t'
                 << (occupancy.empty() ? "1.0" : num_or_null(occupancy)) << '\t'
                 << (u.empty()         ? "\\N" : num_or_null(u))         << '\n';

            PGCopySend(conn, line.str());
        }
    }

    PGCopyEnd(conn);
}

static const char* USAGE_MESSAGE =
    "Usage: cif_parser <input.cif> \"<conninfo>\"\n"
    "Example: parser ccdc-2023-dec.cif "
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

        std::vector<long long> ids = insert_structures(conn, batch, input_cif);
        copy_atoms(conn, batch, ids);

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