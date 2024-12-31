#pragma once

#include <boost/range.hpp>
#include <boost/static_assert.hpp>

namespace fulgor{

typedef uint32_t code_type;

static const size_t max_rules_per_round = 1000;
static const size_t max_dict_size = 1 << 16;
static const size_t min_rule_frequency = 16;

typedef uint64_t int_rule_type;
static const size_t hash_prime = 2013686449;

struct rule_type
{
    rule_type(code_type left, code_type right)
        : m_p((int_rule_type(left) << (8 * sizeof(code_type))) | right)
    {
        // BOOST_STATIC_ASSERT(sizeof(int_rule_type) == 2 * sizeof(code_type));
    }
    
    code_type left() const
    {
        return code_type(m_p >> (8 * sizeof(code_type)));
    }

    code_type right() const
    {
        int_rule_type mask = (int_rule_type(1) << (8 * sizeof(code_type))) - 1;
        return code_type(m_p & mask);
    }

    size_t hash() const
    {
        return m_p * hash_prime;
    }
    
    bool operator==(rule_type const& rhs) const
    {
        return m_p == rhs.m_p;
    }
    
private:
    int_rule_type m_p;
};

const rule_type null_rule(-1,-1); // there can't be 65536 rules

template <typename MappedType>
class rules_table
{
public:
    typedef MappedType mapped_type;
    typedef std::pair<rule_type, mapped_type> value_type;
    typedef typename std::vector<value_type>::iterator iterator;
    
    rules_table()
        : m_table(8, std::make_pair(null_rule, mapped_type()))
    , m_size(0)
    {}

    bool try_get(rule_type const& key, mapped_type& value)
    {
        value_type& cell = get_cell(key);
        if (cell.first == null_rule) {
            return false;
        } else {
            value = cell.second;
            return true;
        }
    }

    mapped_type& operator[](rule_type const& key)
    {
    rehash();
        value_type& cell = get_cell(key);
        
        if (cell.first == null_rule) {
            cell = std::make_pair(key, mapped_type());
    m_size += 1;
        }
        return cell.second;
    }

    iterator begin() { return m_table.begin(); }
    iterator end() { return m_table.end(); }

private:
    
    value_type& get_cell(rule_type const& key)
    {
        size_t h = key.hash();
        size_t mask = m_table.size() - 1;
    value_type* cell = 0;
        while (true) {
    cell = &m_table[h & mask];
            if (cell->first == null_rule ||
                cell->first == key) break;
            ++h;
        }
        return *cell;
    }

    void rehash()
    {
        if (m_table.size() <= m_size * 2) {
            size_t old_size = m_table.size();
    std::vector<value_type> old_table(2 * old_size, std::make_pair(null_rule, mapped_type()));
    old_table.swap(m_table);

            for (size_t i = 0; i < old_table.size(); ++i) {
        if (!(old_table[i].first == null_rule)) {
        get_cell(old_table[i].first) = old_table[i];
        }
            }
        }
    }

    std::vector<value_type> m_table;
size_t m_size;
};

struct second_gt
{
    template <typename Pair>
    bool operator()(Pair const& p1, Pair const& p2) const {
        return p1.second > p2.second;
    }
};

template <typename Range, typename WordVector>
void approximate_repair(Range const& s,
                        std::vector<code_type>& C,
                        WordVector& D,
                        bool preserve_boundaries = false)
{
    typedef typename WordVector::value_type word_type;
    typedef typename word_type::value_type char_type;

    size_t cur_l = boost::size(s);

    // init C
    C.resize(cur_l);
    std::map<char_type, size_t> alph_map;
    alph_map[0] = 1;
    
    typedef typename boost::range_const_iterator<Range>::type iterator_t;
    for (size_t i = 0; i < cur_l; ++i) {
        char_type c = (char_type)*(boost::begin(s) + i);
        size_t code = alph_map[c];
        if (code == 0) {
            code = alph_map.size(); // alph_map already contains c, so this is (old size) + 1
            alph_map[c] = code;
        }
        assert(code - 1 <= std::numeric_limits<code_type>::max());
        C[i] = code_type(code - 1);
    }
    
    // init D
    size_t alph_size = alph_map.size();
    D.resize(alph_size);
    size_t dict_size = D.size();
    for (typename std::map<char_type, size_t>::const_iterator iter = alph_map.begin();
         iter != alph_map.end();
         ++iter) {
        D[iter->second - 1].push_back(iter->first);
        assert(D[iter->second - 1].size() == 1);
    }
 
    std::vector<size_t> L(alph_size, 1);

    // iterate
    typedef rules_table<size_t> map_type;
    map_type counts;
    size_t round = 0;
        
    while (true) {
        for (size_t i = 0; i < cur_l - 1; ++i) {
            // only add rules that would fit in the dictionary
            code_type left = C[i], right = C[i + 1];
            if (dict_size + L[left] + L[right] <= max_dict_size) {
                if (!preserve_boundaries || (left && right)) { // if preserve_boundaries do not create rules that contain 0
                    counts[rule_type(left, right)] += 1;
                }
            }
        }

        // find the best max_rules rules
        std::vector<std::pair<rule_type, size_t> > new_rules;
        for (map_type::iterator iter = counts.begin(); iter != counts.end(); ++iter) {
    if (iter->first == null_rule) continue;
    
            if (iter->second >= min_rule_frequency) {

                if (new_rules.size() < max_rules_per_round) {
                    new_rules.push_back(*iter);
                    std::push_heap(new_rules.begin(), new_rules.end(),
                                   second_gt());
                } else {
                    if (iter->second > new_rules[0].second) {
                        std::pop_heap(new_rules.begin(), new_rules.end(),
                                      second_gt());
                        new_rules.back() = *iter;
                        std::push_heap(new_rules.begin(), new_rules.end(),
                                       second_gt());
                    }
                }
            }
            iter->second = 0; // reset counts
        }

        for (size_t i = new_rules.size(); i > 1; --i) {
            std::pop_heap(new_rules.begin(), new_rules.begin() + i,
                          second_gt());
        }

        if (new_rules.size() == 0) break; // done

        // add to the dictionary all the new rules that fit
        typedef rules_table<code_type> replacements_type;
        replacements_type replacements;

        for (size_t i = 0; i < new_rules.size(); ++i) {
            rule_type const& rule = new_rules[i].first;

            if (dict_size + L[rule.left()] + L[rule.right()] > max_dict_size) {
                continue;
            }

            word_type word(D[rule.left()]);
            word.insert(word.end(), D[rule.right()].begin(), D[rule.right()].end());

            replacements[rule] = code_type(D.size());
            D.push_back(word);
            L.push_back(word.size());
            dict_size += word.size();                    
        }

        // replace all the occurrences in C
        size_t to_i = 0;
        for (size_t from_i = 0; from_i < cur_l;) {
            if (from_i + 2 <= cur_l) {
                rule_type p(C[from_i], C[from_i + 1]);
        code_type new_code;
        if (replacements.try_get(p, new_code)) {
                    C[to_i] = new_code;
                    from_i += 2;
                    to_i += 1;
                    continue;
                }
            }
            // no match found, just copy the current code
            C[to_i++] = C[from_i++];
        }
        cur_l = to_i;
    }
    C.resize(cur_l);
    round += 1;
}

template <typename ColorSets>
struct index<ColorSets>::repair_builder{
    repair_builder() {}

    repair_builder(build_configuration const& build_config): m_build_config(build_config) {}

    void build(index &idx){
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index_already_built");

        index_type index;
        essentials::logger("step 1. loading index to be repair-ed...");
        essentials::load(index, m_build_config.index_filename_to_partition.c_str());
        essentials::logger("DONE");

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        uint64_t num_color_sets = index.num_color_sets();

        {
            essentials::logger("step 2. building repaired colors");
            std::vector<uint32_t> diffs, sizes(num_color_sets);
            typename ColorSets::builder colors_builder(index.num_colors());

            for(uint64_t color_set_id = 0; color_set_id < num_color_sets; ++color_set_id){
                auto it = index.color_set(color_set_id);
                sizes[color_set_id] = it.size();
                
                uint32_t prev = *it;
                ++it;
                diffs.push_back(prev);
                
                for(uint64_t i = 1; i < sizes[color_set_id]; ++i, ++it){
                    diffs.push_back(*it - prev - 1);
                    prev = *it;
                } 
            }

            std::vector<code_type> C;
            std::vector<std::vector<uint32_t>> D;
            approximate_repair(std::make_pair(diffs.data(), diffs.data() + diffs.size()), C, D, false);
            colors_builder.set_code(C);
            colors_builder.set_dict(D);
            colors_builder.build(idx.m_color_sets);
        }  

        {
            essentials::logger("step 3. copy u2c and k2u");
            timer.start();
            idx.m_u2c = index.get_u2c();
            idx.m_k2u = index.get_k2u();
            timer.stop();
            std::cout << "** copying u2c and k2u took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 4. building filenames");
            timer.start();
            idx.m_filenames = index.get_filenames();
            timer.stop();
            std::cout << "** building filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }
    }

private:
    build_configuration m_build_config;
    std::vector<uint64_t> m_offsets;
    std::vector<uint8_t> m_sequences;
};

}
