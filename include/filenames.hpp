#pragma once

#include <cstdint>
#include <fstream>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "external/sshash/external/pthash/external/bits/external/essentials/include/essentials.hpp"

namespace fulgor {

struct filenames {
    void build(std::vector<std::string> const& filenames) {
        uint64_t offset = 0;
        std::vector<uint32_t> offsets;
        std::vector<char> chars;
        offsets.push_back(0);
        for (auto const& f : filenames) {
            std::ranges::copy(f, std::back_inserter(chars));
            offset += f.size();
            offsets.push_back(checked_offset(offset));
        }
        m_offsets = offsets;
        m_chars = chars;
    }

    void build_from_file(std::string const& filepath) {
        std::vector<uint32_t> offsets;
        std::vector<char> chars;
        std::ifstream file(filepath, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open file: " + filepath);
        }
        uint64_t offset = 0;
        offsets.push_back(0);

        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            chars.insert(chars.end(), line.begin(), line.end());

            offset += line.size();
            offsets.push_back(checked_offset(offset));
        }
        m_offsets = offsets;
        m_chars = chars;
    }

    std::string_view operator[](uint64_t i) const {
        uint32_t begin = m_offsets[i];
        uint32_t end = m_offsets[i + 1];
        return {m_chars.data() + begin, end - begin};
    }

    uint64_t num_bits() const {
        return essentials::vec_bytes(m_offsets) * 8 + essentials::vec_bytes(m_chars) * 8;
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

private:
    /* offsets are stored as 32-bit; reject a filename buffer that cannot be addressed */
    static uint32_t checked_offset(const uint64_t offset) {
        if (offset > UINT32_MAX) {
            throw std::runtime_error("reference filenames exceed the 4 GiB offset limit");
        }
        return static_cast<uint32_t>(offset);
    }

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_offsets);
        visitor.visit(t.m_chars);
    }

    essentials::owning_span<uint32_t> m_offsets;
    essentials::owning_span<char> m_chars;
};

}  // namespace fulgor