#include "import_core.hpp"
#include <iostream>
#include <stdexcept>
using namespace unncdb;
void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
int main() {
    try {
        require(sha256("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
                "SHA256");
        auto n = number("1.23(4)e-2");
        require(n.value && std::abs(std::stod(*n.value) - .0123) < 1e-12, "numeric exponent");
        require(n.su && std::abs(std::stod(*n.su) - .0004) < 1e-12, "numeric uncertainty");
        require(number("?").state == "unknown" && number("").state == "absent", "missing states");
        require(!number("NaN").value && !number("Inf").value && !number("1.2garbage").value,
                "invalid numbers");
        std::string source = "# preamble\r\nDATA_a\r\n_tag\r\n;\r\ndata_not_a_block\r\n;\r\n_name "
                             "'data_quoted'\r\n# data_comment\r\n data_b _x 1\n";
        BlockReader r;
        std::vector<RawBlock> blocks;
        std::istringstream in(source);
        r.read(in, [&](RawBlock b) { blocks.push_back(std::move(b)); });
        require(blocks.size() == 2 && blocks[0].name == "a" && blocks[1].name == "b", "framing");
        std::string restored = r.preamble;
        for (auto &b : blocks) {
            require(source.substr(b.offset, b.raw.size()) == b.raw, "offset");
            restored += b.raw;
        }
        require(restored == source, "byte round trip");
        std::string cif =
            "data_test\n_CELL_LENGTH_A 10\n_cell_length_b 10\n_cell_length_c 10\n_cell_angle_alpha "
            "90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n_symmetry_Int_Tables_number "
            "2\nloop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'\nloop_\n_atom_site_label\n_"
            "atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_"
            "atom_site_occupancy\nC1 C 1.2(3)e-1 0 0 .\nC1 C ? 0 0 ?\nCa1? Ca 0 0 0 0.5\n";
        auto p = parse_block({1, 0, "test", cif});
        require(p.atoms.size() == 3, "keep all atom rows");
        require(p.atoms[0]["occupancy_state"] == Value{"not_applicable"}, "dot state");
        require(!p.atoms[1]["x"] && p.fields["extraction_status"] == Value{"issues"},
                "missing coords flagged");
        require(p.fields["is_centrosymmetric"] == Value{"true"}, "centrosymmetry");
        require(p.fields["space_group_num"] == Value{"2"}, "space group");
        bool duplicate = false, ambiguous = false;
        for (auto &i : p.issues) {
            duplicate |= i.code == "duplicate_atom_label";
            ambiguous |= i.code == "ambiguous_atom_label";
        }
        require(duplicate && ambiguous, "quality diagnostics");
        auto bad = parse_block({1, 0, "bad", "data_bad\n_tag\n"});
        require(bad.parse_status == "parse_error", "malformed block");
        auto conflict_cif = cif;
        auto declared_pos = conflict_cif.find("_symmetry_Int_Tables_number 2");
        conflict_cif.replace(declared_pos, std::string("_symmetry_Int_Tables_number 2").size(),
                             "_symmetry_Int_Tables_number 1");
        auto conflict = parse_block({1, 0, "conflict", conflict_cif});
        require(!conflict.fields["is_centrosymmetric"] &&
                    conflict.fields["symmetry_status"] == Value{"conflict"},
                "conflicting labels withheld");
        auto missing_ops = cif;
        auto ops_pos = missing_ops.find("loop_\n_symmetry_equiv_pos_as_xyz");
        auto ops_end = missing_ops.find("loop_\n_atom_site_label", ops_pos);
        missing_ops.erase(ops_pos, ops_end - ops_pos);
        auto declared_only = parse_block({1, 0, "declared", missing_ops});
        require(declared_only.fields["is_centrosymmetric"] == Value{"true"} &&
                    declared_only.fields["symmetry_status"] == Value{"declared_only"},
                "declared-only symmetry");
        auto invalid_utf8 =
            parse_block({1, 0, "bad_utf8", std::string("data_bad\n_tag ") + char(0xff) + "\n"});
        require(invalid_utf8.parse_status == "parse_error", "invalid UTF8 diagnosed");
        std::istringstream inline_blocks("data_one _tag 'data_quoted' data_two _tag 2");
        BlockReader inline_reader;
        size_t inline_count = 0;
        inline_reader.read(inline_blocks, [&](RawBlock) { ++inline_count; });
        require(inline_count == 2, "inline headers and no final newline");
        require(element_symbol({}, Value{"Ca1"}) == Value{"Ca"}, "element fallback");
        require(!element_symbol(Value{"Bogus"}, Value{"C1"}),
                "invalid supplied element not guessed");
        std::cout << "All import core checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
