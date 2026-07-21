#pragma once

namespace fulgor {

struct filenames {
    void build(std::vector<std::string> const& filenames) {
        uint32_t offset = 0;
        std::vector<uint32_t> offsets;
        std::vector<char> chars;
        offsets.push_back(offset);
        for (auto const& f : filenames) {
            std::ranges::copy(f, std::back_inserter(chars));
            offset += f.size();
            offsets.push_back(offset);
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
        uint32_t offset = 0;
        offsets.push_back(offset);

        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            chars.insert(chars.end(), line.begin(), line.end());

            offset += line.size();
            offsets.push_back(offset);
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
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_offsets);
        visitor.visit(t.m_chars);
    }

    essentials::owning_span<uint32_t> m_offsets;
    essentials::owning_span<char> m_chars;
};

}  // namespace fulgor