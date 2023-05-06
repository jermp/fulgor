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

    std::string_view filename(uint64_t doc_id) const {
        uint32_t begin = m_offsets[doc_id];
        uint32_t end = m_offsets[doc_id + 1];
        return {m_chars.data() + begin, end - begin};
    }

    // void print() const {
    //     for (uint64_t i = 0; i != num_docs(); ++i) {
    //         auto const& s = filename(i);
    //         std::cout << s << std::endl;
    //     }
    // }

    uint32_t num_docs() const {
        assert(m_offsets.size() > 0);
        return m_offsets.size() - 1;
    }

    uint64_t num_bits() const {
        return essentials::vec_bytes(m_offsets) * 8 + essentials::vec_bytes(m_chars) * 8;
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_offsets);
        visitor.visit(m_chars);
    }

private:
    std::vector<uint32_t> m_offsets;
    std::vector<char> m_chars;
};

}  // namespace fulgor