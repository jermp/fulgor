#pragma once

namespace fulgor {

struct filenames {
    void build(std::vector<std::string> const& filenames) {
        uint32_t offset = 0;
        m_offsets.push_back(offset);
        for (auto const& f : filenames) {
            std::copy(f.begin(), f.end(), std::back_inserter(m_chars));
            offset += f.size();
            m_offsets.push_back(offset);
        }
    }

    std::string_view operator[](uint64_t i) const {
        uint32_t begin = m_offsets[i];
        uint32_t end = m_offsets[i + 1];
        return {m_chars.data() + begin, end - begin};
    }

    // uint32_t num_docs() const {
    //     assert(m_offsets.size() > 0);
    //     return m_offsets.size() - 1;
    // }

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

    std::vector<uint32_t> m_offsets;
    std::vector<char> m_chars;
};

}  // namespace fulgor