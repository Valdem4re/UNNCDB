#pragma once
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <fstream>
#include <functional>
#include <gemmi/cif.hpp>
#include <gemmi/elem.hpp>
#include <gemmi/numb.hpp>
#include <gemmi/symmetry.hpp>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <openssl/evp.h>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace unncdb {
inline constexpr const char *VERSION = "2.0.1";
inline bool compatible_version(const std::string &version) {
    // 2.0.1 changes performance only; derived database values are unchanged.
    return version == "2.0.0" || version == VERSION;
}
using Value = std::optional<std::string>;
using Row = std::map<std::string, Value>;
inline std::string lower(std::string s) {
    for (char &c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}
inline std::string decimal(double v) {
    std::array<char, 64> buffer;
    auto result =
        std::to_chars(buffer.data(), buffer.data() + buffer.size(), v, std::chars_format::general,
                      std::numeric_limits<double>::max_digits10);
    if (result.ec != std::errc())
        throw std::runtime_error("Cannot format numeric value");
    return std::string(buffer.data(), result.ptr);
}

class Sha256 {
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> ctx{EVP_MD_CTX_new(), EVP_MD_CTX_free};

  public:
    Sha256() {
        if (!ctx || EVP_DigestInit_ex(ctx.get(), EVP_sha256(), nullptr) != 1)
            throw std::runtime_error("SHA256 init failed");
    }
    void update(const std::string &data) {
        if (EVP_DigestUpdate(ctx.get(), data.data(), data.size()) != 1)
            throw std::runtime_error("SHA256 update failed");
    }
    std::string finish() {
        unsigned char bytes[EVP_MAX_MD_SIZE];
        unsigned int n = 0;
        if (EVP_DigestFinal_ex(ctx.get(), bytes, &n) != 1)
            throw std::runtime_error("SHA256 finish failed");
        std::ostringstream s;
        s << std::hex << std::setfill('0');
        for (unsigned int i = 0; i < n; ++i)
            s << std::setw(2) << static_cast<unsigned>(bytes[i]);
        return s.str();
    }
};
inline std::string sha256(const std::string &s) {
    Sha256 hash;
    hash.update(s);
    return hash.finish();
}
inline std::string file_sha256(const std::string &path) {
    std::ifstream f(path, std::ios::binary);
    if (!f)
        throw std::runtime_error("Cannot open input: " + path);
    Sha256 hash;
    std::array<char, 1024 * 1024> buf;
    while (f) {
        f.read(buf.data(), buf.size());
        if (f.gcount())
            hash.update(std::string(buf.data(), static_cast<size_t>(f.gcount())));
    }
    if (!f.eof())
        throw std::runtime_error("Cannot read input: " + path);
    return hash.finish();
}
struct RawBlock {
    long long ordinal = 0, offset = 0;
    std::string name, raw;
};
// Framing only: Gemmi performs syntax validation. Semicolon text fields,
// quoted tokens and comments cannot introduce data-block boundaries.
// Preamble + concatenated raw blocks reconstruct the original bytes.
class BlockReader {
    bool text_field = false;
    static bool ws(char c) {
        return c == ' ' || c == '\t' || c == '\r' || c == '\n';
    }
    std::vector<std::pair<size_t, std::string>> headers(const std::string &line) {
        std::vector<std::pair<size_t, std::string>> result;
        if (!line.empty() && line[0] == ';') {
            text_field = !text_field;
            return result;
        }
        if (text_field)
            return result;
        size_t p = 0;
        while (p < line.size()) {
            while (p < line.size() && ws(line[p]))
                ++p;
            if (p == line.size() || line[p] == '#')
                break;
            if (line[p] == '\'' || line[p] == '"') {
                char quote = line[p++];
                while (p < line.size()) {
                    if (line[p] == quote && (p + 1 == line.size() || ws(line[p + 1]))) {
                        ++p;
                        break;
                    }
                    ++p;
                }
            } else {
                size_t start = p;
                while (p < line.size() && !ws(line[p]))
                    ++p;
                std::string word = line.substr(start, p - start);
                if (word.size() > 5 && lower(word.substr(0, 5)) == "data_")
                    result.emplace_back(start, word.substr(5));
            }
        }
        return result;
    }

  public:
    std::string preamble;
    template <class Callback> void read(std::istream &in, Callback callback) {
        RawBlock current;
        std::string line;
        long long offset = 0, ordinal = 0;
        while (std::getline(in, line)) {
            if (!in.eof())
                line.push_back('\n');
            auto found = headers(line);
            size_t start = 0;
            for (const auto &[pos, name] : found) {
                (current.ordinal ? current.raw : preamble).append(line, start, pos - start);
                if (current.ordinal)
                    callback(std::move(current));
                current = RawBlock{++ordinal, offset + static_cast<long long>(pos), name, {}};
                start = pos;
            }
            (current.ordinal ? current.raw : preamble).append(line, start, std::string::npos);
            offset += static_cast<long long>(line.size());
        }
        if (!in.eof())
            throw std::runtime_error("Input read failed");
        if (current.ordinal)
            callback(std::move(current));
    }
};
struct Number {
    Value value, su;
    std::string state;
};
inline Number number(const std::string &raw) {
    if (raw.empty())
        return {{}, {}, "absent"};
    if (raw == "?")
        return {{}, {}, "unknown"};
    if (raw == ".")
        return {{}, {}, "not_applicable"};
    if (raw.find('(') == std::string::npos) {
        double v = gemmi::cif::as_number(raw);
        return std::isfinite(v) ? Number{decimal(v), {}, "explicit"} : Number{{}, {}, "invalid"};
    }
    // Preserve exponent even when s.u. precedes it; accept s.u. after exponent too.
    static const std::regex pattern(
        R"(^([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+))(?:\(([0-9]+)\))?([eE][+-]?[0-9]+)?(?:\(([0-9]+)\))?$)");
    std::smatch m;
    if (!std::regex_match(raw, m, pattern) || (m[2].matched && m[4].matched))
        return {{}, {}, "invalid"};
    try {
        const std::string mantissa = m[1], exponent = m[3];
        double value = std::stod(mantissa + exponent);
        if (!std::isfinite(value))
            return {{}, {}, "invalid"};
        Value su;
        std::string digits = m[2].matched ? m[2].str() : m[4].str();
        if (!digits.empty()) {
            auto dot = mantissa.find('.');
            int places = dot == std::string::npos ? 0 : static_cast<int>(mantissa.size() - dot - 1);
            int exp = exponent.empty() ? 0 : std::stoi(exponent.substr(1));
            double u = std::stod(digits) * std::pow(10., exp - places);
            if (!std::isfinite(u))
                return {{}, {}, "invalid"};
            su = decimal(u);
        }
        return {decimal(value), su, "explicit"};
    } catch (...) {
        return {{}, {}, "invalid"};
    }
}
inline Value text(const std::string &raw) {
    if (raw.empty() || gemmi::cif::is_null(raw))
        return {};
    return gemmi::cif::as_string(raw);
}
struct Issue {
    std::string severity, code;
    Value row, tag;
    std::string message;
};
struct Parsed {
    std::string parse_status = "parsed";
    bool structure = false;
    Row fields;
    std::vector<Row> atoms, adp, operations;
    std::vector<Issue> issues;
    void issue(std::string code, std::string msg, Value row = {}, Value tag = {},
               std::string severity = "error") {
        issues.push_back(
            {std::move(severity), std::move(code), std::move(row), std::move(tag), std::move(msg)});
    }
};
inline std::string scalar(gemmi::cif::Block &b, std::initializer_list<const char *> tags) {
    for (auto tag : tags)
        if (auto *p = b.find_value(tag))
            if (!gemmi::cif::is_null(*p))
                return *p;
    return {};
}
// Column lookup is linear in CIF tags. Resolve once per block, not once per atom.
class ColumnCache {
    gemmi::cif::Block &block;
    std::map<std::string, gemmi::cif::Column, std::less<>> columns;

  public:
    explicit ColumnCache(gemmi::cif::Block &b) : block(b) {}
    std::string get(const char *tag, int row) {
        auto it = columns.find(tag);
        if (it == columns.end())
            it = columns.emplace(tag, block.find_values(tag)).first;
        const auto &c = it->second;
        return row < static_cast<int>(c.length()) ? c.at(row) : "";
    }
};
// Cache complete, canonical operation sets, never just the declared group number.
inline const gemmi::SpaceGroup *cached_spacegroup(const gemmi::GroupOps &group) {
    thread_local std::map<std::vector<gemmi::Op>, const gemmi::SpaceGroup *> cache;
    auto key = group.all_ops_sorted();
    auto found = cache.find(key);
    if (found != cache.end())
        return found->second;
    auto *result = gemmi::find_spacegroup_by_ops(group);
    if (cache.size() < 4096)
        cache.emplace(std::move(key), result);
    return result;
}
inline Value numeric(Parsed &p, const std::string &raw, const std::string &tag, Value row = {}) {
    auto n = number(raw);
    if (n.state == "invalid")
        p.issue("invalid_number", "Invalid numeric token: " + raw, row, tag);
    return n.value;
}
inline std::optional<int> integer(const Value &n) {
    if (!n)
        return {};
    double v = std::stod(*n);
    if (v != std::floor(v) || v < INT32_MIN || v > INT32_MAX)
        return {};
    return static_cast<int>(v);
}
inline Value element_symbol(const Value &type, const Value &label) {
    auto parse = [](const std::string &s, bool from_label) -> Value {
        if (s.empty() || !std::isupper(static_cast<unsigned char>(s[0])))
            return {};
        size_t n = 1;
        if (s.size() > 1 && std::islower(static_cast<unsigned char>(s[1])))
            n = 2;
        auto symbol = s.substr(0, n);
        if (symbol == "T")
            return symbol;
        if (gemmi::Element(symbol).elem == gemmi::El::X)
            return {};
        if (!from_label)
            for (size_t i = n; i < s.size(); ++i)
                if (!std::isdigit(static_cast<unsigned char>(s[i])) && s[i] != '+' && s[i] != '-')
                    return {};
        return symbol;
    };
    // A supplied but invalid type must not be silently replaced by a label guess.
    return type ? parse(*type, false) : (label ? parse(*label, true) : Value{});
}
inline void extract_symmetry(gemmi::cif::Block &b, Parsed &p) {
    auto raw_num = scalar(b, {"_space_group_IT_number", "_symmetry_Int_Tables_number"});
    p.fields["space_group_num_raw"] = raw_num.empty() ? Value{} : Value{raw_num};
    auto declared = integer(numeric(p, raw_num, "_space_group_IT_number"));
    if (declared && (*declared < 1 || *declared > 230)) {
        p.issue("invalid_space_group_number", "Declared space-group number is outside 1..230");
        declared.reset();
    } else if (!raw_num.empty() && !declared)
        p.issue("invalid_space_group_number", "Declared space-group number is not an integer");
    if (declared)
        p.fields["declared_space_group_num"] = std::to_string(*declared);
    auto hm = text(scalar(b, {"_space_group_name_H-M_alt", "_symmetry_space_group_name_H-M"}));
    auto hall = text(scalar(b, {"_space_group_name_Hall", "_symmetry_space_group_name_Hall"}));
    p.fields["space_group"] = hm;
    p.fields["hall_symbol"] = hall;
    const gemmi::SpaceGroup *by_name = nullptr, *by_ops = nullptr;
    bool conflict = false;
    std::vector<gemmi::Op> ops;
    const char *op_tag = b.has_tag("_space_group_symop_operation_xyz")
                             ? "_space_group_symop_operation_xyz"
                             : "_symmetry_equiv_pos_as_xyz";
    auto column = b.find_values(op_tag);
    for (int i = 0; i < static_cast<int>(column.length()); ++i) {
        auto xyz = text(column.at(i));
        p.operations.push_back({{"source_row", std::to_string(i + 1)},
                                {"operation_xyz", xyz ? xyz : Value{column.at(i)}}});
        try {
            if (!xyz)
                throw std::runtime_error("Unknown operation");
            ops.push_back(gemmi::parse_triplet(*xyz).wrap());
        } catch (const std::exception &e) {
            p.issue("invalid_symmetry_operation", e.what(), std::to_string(i + 1), op_tag);
            conflict = true;
        }
    }
    try {
        if (hm)
            by_name = gemmi::find_spacegroup_by_name(
                *hm, p.fields["alpha"] ? std::stod(*p.fields["alpha"]) : 0.,
                p.fields["gamma"] ? std::stod(*p.fields["gamma"]) : 0.);
        if (hm && !by_name)
            p.issue("unknown_space_group_name", "Cannot resolve H-M symbol", {}, {}, "warning");
        if (!ops.empty() && !conflict) {
            std::set<gemmi::Op> unique(ops.begin(), ops.end());
            if (unique.size() != ops.size())
                p.issue("duplicate_symmetry_operation", "Repeated symmetry operation", {}, {},
                        "warning");
            std::vector<gemmi::Op> distinct(unique.begin(), unique.end());
            auto group = gemmi::split_centering_vectors(distinct);
            by_ops = cached_spacegroup(group);
            if (!by_ops) {
                p.issue("unresolved_symmetry_operations",
                        "Operations do not match a supported complete group in this setting");
                conflict = true;
            }
        }
        if (hall) {
            auto hall_ops = gemmi::symops_from_hall(hall->c_str());
            auto by_hall = cached_spacegroup(hall_ops);
            if (!by_hall) {
                p.issue("unresolved_hall_symbol", "Hall symbol could not be matched");
                conflict = true;
            }
            if (by_ops && !hall_ops.is_same_as(by_ops->operations())) {
                p.issue("symmetry_conflict", "Hall symbol disagrees with supplied operations");
                conflict = true;
            }
            if (by_name && by_hall && by_name->number != by_hall->number) {
                p.issue("symmetry_conflict", "H-M and Hall symbols disagree");
                conflict = true;
            }
            if (!by_ops && column.length() == 0)
                by_ops = by_hall;
        }
        const gemmi::SpaceGroup *resolved = by_ops ? by_ops : by_name;
        if (resolved && declared && resolved->number != *declared) {
            p.issue("symmetry_conflict", "Declared number disagrees with resolved group");
            conflict = true;
        }
        if (by_ops && by_name && by_ops->number != by_name->number) {
            p.issue("symmetry_conflict", "H-M symbol disagrees with operations");
            conflict = true;
        }
        if (!resolved && declared)
            resolved = &gemmi::get_spacegroup_reference_setting(*declared);
        if (resolved && !conflict) {
            p.fields["space_group_num"] = std::to_string(resolved->number);
            p.fields["point_group"] = resolved->point_group_hm();
            p.fields["crystal_system"] = resolved->crystal_system_str();
            p.fields["is_centrosymmetric"] = resolved->is_centrosymmetric() ? "true" : "false";
            p.fields["symmetry_status"] = by_ops ? "verified" : "declared_only";
            p.fields["symmetry_method"] = column.length() ? "operations"
                                          : hall          ? "hall"
                                          : by_name       ? "hm"
                                                          : "number";
        } else {
            p.fields["symmetry_status"] = conflict ? "conflict" : "unknown";
            if (!resolved && !conflict)
                p.issue("missing_symmetry", "Cannot determine space-group type");
        }
        p.fields["symmetry_version"] = "gemmi-0.7.5";
    } catch (const std::exception &e) {
        p.fields["symmetry_status"] = "conflict";
        p.issue("symmetry_error", e.what());
    }
}
inline bool valid_utf8_text(const std::string &s) {
    size_t i = 0;
    while (i < s.size()) {
        auto c = static_cast<unsigned char>(s[i++]);
        if (c == 0)
            return false;
        if (c < 0x80)
            continue;
        int n;
        unsigned code, min;
        if (c >= 0xc2 && c <= 0xdf) {
            n = 1;
            code = c & 31;
            min = 0x80;
        } else if (c >= 0xe0 && c <= 0xef) {
            n = 2;
            code = c & 15;
            min = 0x800;
        } else if (c >= 0xf0 && c <= 0xf4) {
            n = 3;
            code = c & 7;
            min = 0x10000;
        } else
            return false;
        if (i + static_cast<size_t>(n) > s.size())
            return false;
        while (n--) {
            auto d = static_cast<unsigned char>(s[i++]);
            if ((d & 0xc0) != 0x80)
                return false;
            code = (code << 6) | (d & 63);
        }
        if (code < min || code > 0x10ffff || (code >= 0xd800 && code <= 0xdfff))
            return false;
    }
    return true;
}
inline Parsed parse_block(const RawBlock &raw) {
    Parsed p;
    if (!valid_utf8_text(raw.raw)) {
        p.parse_status = "parse_error";
        p.issue("invalid_text_encoding", "Block is not NUL-free UTF-8; exact bytes retained");
        return p;
    }
    try {
        auto doc = gemmi::cif::read_memory(raw.raw.data(), raw.raw.size(), raw.name.c_str(), 2);
        if (doc.blocks.size() != 1)
            throw std::runtime_error("Expected exactly one CIF data block");
        auto &b = doc.blocks.front();
        auto labels = b.find_values("_atom_site_label");
        auto x = b.find_values("_atom_site_fract_x"), y = b.find_values("_atom_site_fract_y"),
             z = b.find_values("_atom_site_fract_z");
        size_t count = std::max({labels.length(), x.length(), y.length(), z.length(),
                                 b.find_values("_atom_site_type_symbol").length()});
        p.structure = count != 0;
        if (!p.structure)
            return p;
        // Atomic columns must be in the same source table; do not zip unrelated loops.
        gemmi::cif::Item *atom_table = nullptr;
        for (const char *tag :
             {"_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
              "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy"}) {
            auto c = b.find_values(tag);
            if (c.length()) {
                if (!atom_table)
                    atom_table = c.item();
                else if (c.item() != atom_table) {
                    p.issue("split_atom_table", "Atom-site columns are not in one table", {}, tag);
                }
            }
        }
        p.fields["refcode"] = text(scalar(b, {"_database_code_CSD"}));
        if (!p.fields["refcode"])
            p.fields["refcode"] = b.name;
        p.fields["formula"] = text(scalar(b, {"_chemical_formula_sum"}));
        p.fields["formula_moiety"] = text(scalar(b, {"_chemical_formula_moiety"}));
        p.fields["description"] =
            text(scalar(b, {"_chemical_name_systematic", "_chemical_name_common"}));
        p.fields["refinement_details"] = text(scalar(b, {"_refine_special_details"}));
        const std::vector<std::pair<const char *, const char *>> numbers = {
            {"a", "_cell_length_a"},
            {"b", "_cell_length_b"},
            {"c", "_cell_length_c"},
            {"alpha", "_cell_angle_alpha"},
            {"beta", "_cell_angle_beta"},
            {"gamma", "_cell_angle_gamma"},
            {"volume", "_cell_volume"},
            {"r_factor_gt", "_refine_ls_R_factor_gt"},
            {"r_factor_all", "_refine_ls_R_factor_all"},
            {"wr_factor_gt", "_refine_ls_wR_factor_gt"},
            {"wr_factor_ref", "_refine_ls_wR_factor_ref"}};
        for (auto [field, tag] : numbers) {
            auto raw_number = scalar(b, {tag});
            p.fields[field] = numeric(p, raw_number, tag);
            if (std::string(field) == "a" || std::string(field) == "b" ||
                std::string(field) == "c" || std::string(field) == "alpha" ||
                std::string(field) == "beta" || std::string(field) == "gamma")
                p.fields[std::string(field) + "_su"] = number(raw_number).su;
        }
        p.fields["radiation_probe"] = text(scalar(b, {"_diffrn_radiation_probe"}));
        auto wavelengths = b.find_values("_diffrn_radiation_wavelength");
        if (wavelengths.length() == 1)
            p.fields["radiation_wavelength"] =
                numeric(p, wavelengths.at(0), "_diffrn_radiation_wavelength");
        else if (wavelengths.length() > 1)
            p.issue("multiple_wavelengths", "Multiple wavelengths retained in raw CIF", {}, {},
                    "warning");
        p.fields["measurement_temperature"] =
            numeric(p, scalar(b, {"_diffrn_ambient_temperature", "_cell_measurement_temperature"}),
                    "_diffrn_ambient_temperature");
        auto zn =
            integer(numeric(p, scalar(b, {"_cell_formula_units_Z"}), "_cell_formula_units_Z"));
        if (zn && *zn > 0) {
            p.fields["z"] = std::to_string(*zn);
            if (p.fields["volume"])
                p.fields["reduced_volume"] = decimal(std::stod(*p.fields["volume"]) / *zn);
        }
        bool valid_cell = true;
        for (const char *key : {"a", "b", "c", "alpha", "beta", "gamma"}) {
            auto v = p.fields[key];
            if (!v || std::stod(*v) <= 0 || (std::string(key).size() > 1 && std::stod(*v) >= 180))
                valid_cell = false;
        }
        if (valid_cell) {
            double ca = std::cos(std::stod(*p.fields["alpha"]) * std::acos(-1.) / 180.),
                   cb = std::cos(std::stod(*p.fields["beta"]) * std::acos(-1.) / 180.),
                   cg = std::cos(std::stod(*p.fields["gamma"]) * std::acos(-1.) / 180.);
            valid_cell = 1 + 2 * ca * cb * cg - ca * ca - cb * cb - cg * cg > 1e-12;
        }
        if (!valid_cell)
            p.issue("invalid_cell", "Missing or nonphysical cell parameters");
        if (p.fields["r_factor_gt"] &&
            (std::stod(*p.fields["r_factor_gt"]) < 0 || std::stod(*p.fields["r_factor_gt"]) > 1))
            p.issue("unusual_r_factor", "R-factor is outside [0,1]; original value retained");
        extract_symmetry(b, p);
        bool disorder = false;
        for (const char *key : {"description", "refinement_details"})
            if (p.fields[key]) {
                static const std::regex keyword(R"(\bdisorder(?:ed|ing)?\b)", std::regex::icase);
                if (std::regex_search(*p.fields[key], keyword)) {
                    disorder = true;
                    p.issue("disorder_text", "Description mentions disorder; inspect context", {},
                            {}, "warning");
                }
            }
        ColumnCache cells(b);
        std::map<std::string, int> label_count;
        p.atoms.reserve(count);
        for (size_t i = 0; i < count; ++i) {
            Row a;
            Value row = std::to_string(i + 1);
            a["source_row"] = row;
            a["label"] = text(cells.get("_atom_site_label", i));
            a["raw_type_symbol"] = text(cells.get("_atom_site_type_symbol", i));
            a["element"] = element_symbol(a["raw_type_symbol"], a["label"]);
            a["element_source"] = a["raw_type_symbol"] ? "type_symbol" : "label";
            if (!a["label"])
                p.issue("missing_atom_label", "Missing atom label", row);
            else {
                if (++label_count[*a["label"]] > 1)
                    p.issue("duplicate_atom_label", "Repeated label: " + *a["label"], row);
                if (a["label"]->find('?') != std::string::npos)
                    p.issue("ambiguous_atom_label", "Question mark in atom label", row);
            }
            if (!a["element"])
                p.issue("unknown_element", "Cannot determine atom element", row);
            for (const char *axis : {"x", "y", "z"}) {
                std::string tag = std::string("_atom_site_fract_") + axis;
                auto n = number(cells.get(tag.c_str(), i));
                a[axis] = n.value;
                a[std::string(axis) + "_su"] = n.su;
                if (!n.value)
                    p.issue("invalid_atom_coordinate", "Missing or invalid fractional coordinate",
                            row, tag);
            }
            auto occ_raw = cells.get("_atom_site_occupancy", i);
            auto occ = number(occ_raw);
            a["occupancy_raw"] = occ_raw.empty() ? Value{} : Value{occ_raw};
            a["occupancy_state"] = occ.state;
            a["occupancy"] = occ.state == "absent" ? Value{"1"} : occ.value;
            if (occ.state == "invalid")
                p.issue("invalid_occupancy", "Invalid occupancy", row);
            if (occ.value && (std::stod(*occ.value) < 0 || std::stod(*occ.value) > 1))
                p.issue("unusual_occupancy", "Occupancy outside [0,1]", row);
            if (occ.value && std::stod(*occ.value) > 0 && std::stod(*occ.value) < 1)
                p.issue("partial_occupancy",
                        "Partial occupancy; may indicate disorder or non-stoichiometry", row, {},
                        "warning");
            for (auto [field, tag] : std::vector<std::pair<const char *, const char *>>{
                     {"u_iso", "_atom_site_U_iso_or_equiv"},
                     {"b_iso", "_atom_site_B_iso_or_equiv"}})
                a[field] = numeric(p, cells.get(tag, i), tag, row);
            for (auto [field, tag] : std::vector<std::pair<const char *, const char *>>{
                     {"adp_type", "_atom_site_adp_type"},
                     {"disorder_assembly", "_atom_site_disorder_assembly"},
                     {"disorder_group", "_atom_site_disorder_group"},
                     {"refinement_flags", "_atom_site_refinement_flags"},
                     {"calc_flag", "_atom_site_calc_flag"}})
                a[field] = text(cells.get(tag, i));
            if (a["disorder_group"] && *a["disorder_group"] != "0")
                disorder = true;
            if (a["disorder_assembly"] && *a["disorder_assembly"] != "0")
                disorder = true;
            p.atoms.push_back(std::move(a));
        }
        // Keep every ADP row, even when labels cannot be joined unambiguously.
        auto aniso = b.find_values("_atom_site_aniso_label");
        for (int i = 0; i < static_cast<int>(aniso.length()); ++i) {
            Row a;
            auto row = std::to_string(i + 1);
            a["source_row"] = row;
            a["label"] = text(aniso.at(i));
            if (!a["label"] || label_count[*a["label"]] != 1)
                p.issue("ambiguous_adp_link", "ADP label has zero or multiple matching sites", row);
            bool has_u = false, has_b = false;
            for (const char *term : {"11", "22", "33", "12", "13", "23"})
                for (const char *form : {"U", "B"}) {
                    std::string tag = std::string("_atom_site_aniso_") + form + "_" + term;
                    auto raw_value = cells.get(tag.c_str(), i);
                    auto field = lower(form) + term;
                    a[field] = numeric(p, raw_value, tag, row);
                    if (!raw_value.empty()) {
                        if (*form == 'U')
                            has_u = true;
                        else
                            has_b = true;
                    }
                }
            a["representation"] = has_u ? (has_b ? "U+B" : "U") : (has_b ? "B" : "unknown");
            p.adp.push_back(std::move(a));
        }
        p.fields["disorder_status"] = disorder ? "suspected" : "unknown";
        p.fields["atom_rows_source"] = std::to_string(count);
        p.fields["atom_rows_stored"] = std::to_string(p.atoms.size());
        bool errors = false;
        for (const auto &issue : p.issues)
            if (issue.severity == "error")
                errors = true;
        p.fields["extraction_status"] = errors ? "issues" : "ok";
    } catch (const std::exception &e) {
        p.parse_status = p.structure ? "extraction_error" : "parse_error";
        p.issue(p.structure ? "extraction_error" : "parse_error", e.what());
        if (p.structure)
            p.fields["extraction_status"] = "issues";
    }
    return p;
}
} // namespace unncdb
