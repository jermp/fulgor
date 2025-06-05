#include "include/color_sets/hybrid.hpp"
#include "include/color_sets/differential.hpp"
#include "include/color_sets/meta.hpp"
#include "include/color_sets/meta_differential.hpp"

namespace fulgor {

void hybrid::print_stats() const  //
{
    const uint64_t num_buckets = 10;
    assert(num_buckets > 0);
    uint64_t bucket_size = m_num_colors / num_buckets;
    std::vector<uint32_t> color_set_size_upperbounds;
    for (uint64_t i = 0, curr_color_set_size_upper_bound = bucket_size; i != num_buckets;
         ++i, curr_color_set_size_upper_bound += bucket_size) {
        if (i == num_buckets - 1) curr_color_set_size_upper_bound = m_num_colors;
        color_set_size_upperbounds.push_back(curr_color_set_size_upper_bound);
    }

    std::vector<uint64_t> num_bits_per_bucket;
    std::vector<uint64_t> num_color_sets_per_bucket;
    std::vector<uint64_t> num_ints_per_bucket;
    num_bits_per_bucket.resize(num_buckets, 0);
    num_color_sets_per_bucket.resize(num_buckets, 0);
    num_ints_per_bucket.resize(num_buckets, 0);

    uint64_t num_total_integers = 0;
    for (uint64_t color_set_id = 0; color_set_id != m_offsets.size() - 1; ++color_set_id) {
        uint64_t offset = m_offsets.access(color_set_id);
        auto it = m_color_sets.get_iterator_at(offset);
        uint32_t color_set_size = bits::util::read_delta(it);
        uint64_t num_bits = m_offsets.access(color_set_id + 1) - offset;
        auto bucket_it = std::upper_bound(color_set_size_upperbounds.begin(),
                                          color_set_size_upperbounds.end(), color_set_size);
        if (bucket_it != color_set_size_upperbounds.begin() and
            *(bucket_it - 1) == color_set_size) {
            --bucket_it;
        }
        uint64_t bucket_index = std::distance(color_set_size_upperbounds.begin(), bucket_it);
        num_bits_per_bucket[bucket_index] += num_bits;
        num_color_sets_per_bucket[bucket_index] += 1;
        num_ints_per_bucket[bucket_index] += color_set_size;
        num_total_integers += color_set_size;
    }

    std::cout << "Color sets space breakdown:\n";
    uint64_t integers = 0;
    uint64_t bits = 0;
    const uint64_t total_bits = num_bits();
    for (uint64_t i = 0, curr_color_set_size_upper_bound = 0; i != num_buckets; ++i) {
        if (i == num_buckets - 1) {
            curr_color_set_size_upper_bound = m_num_colors;
        } else {
            curr_color_set_size_upper_bound += bucket_size;
        }
        if (num_color_sets_per_bucket[i] > 0) {
            uint64_t n = num_ints_per_bucket[i];
            integers += n;
            bits += num_bits_per_bucket[i];
            std::cout << "  num. color_sets of size > "
                      << (curr_color_set_size_upper_bound - bucket_size)
                      << " and <= " << curr_color_set_size_upper_bound << ": "
                      << num_color_sets_per_bucket[i] << " ("
                      << (num_color_sets_per_bucket[i] * 100.0) / num_color_sets()
                      << "%) -- integers: " << n << " (" << (n * 100.0) / num_total_integers
                      << "%) -- bits/int: " << static_cast<double>(num_bits_per_bucket[i]) / n
                      << " -- " << static_cast<double>(num_bits_per_bucket[i]) / total_bits * 100.0
                      << "\% of total space" << '\n';
        }
    }
    assert(integers == num_total_integers);
    assert(std::accumulate(num_color_sets_per_bucket.begin(), num_color_sets_per_bucket.end(),
                           uint64_t(0)) == num_color_sets());
    std::cout << "  colors: " << static_cast<double>(bits) / integers << " bits/int" << std::endl;
    std::cout << "  offsets: "
              << ((sizeof(m_num_colors) + sizeof(m_sparse_set_threshold_size) +
                   sizeof(m_very_dense_set_threshold_size) + m_offsets.num_bytes()) *
                  8.0) /
                     integers
              << " bits/int" << std::endl;
}

template <typename ColorSets>
void meta<ColorSets>::print_stats() const  //
{
    std::cout << "Color sets statistics:\n";
    std::cout << "  Number of partitions: " << num_partitions() << '\n';
    uint64_t num_bits_colors = 0;

    uint64_t num_partial_color_sets_very_dense = 0;
    uint64_t num_partial_color_sets_dense = 0;
    uint64_t num_partial_color_sets_sparse = 0;
    uint64_t num_total_partial_colors = 0;

    for (auto const& pcs : m_partial_color_sets) {
        // pcs.print_stats();
        const uint64_t n = pcs.num_color_sets();
        num_total_partial_colors += n;
        for (uint64_t i = 0; i != n; ++i) {
            auto it = pcs.color_set(i);
            if (it.encoding_type() == encoding_t::complement_delta_gaps) {
                ++num_partial_color_sets_very_dense;
            } else if (it.encoding_type() == encoding_t::bitmap) {
                ++num_partial_color_sets_dense;
            } else {
                assert(it.encoding_type() == encoding_t::delta_gaps);
                ++num_partial_color_sets_sparse;
            }
        }

        num_bits_colors += pcs.num_bits();
    }

    assert(num_total_partial_colors > 0);
    assert(num_bits() > 0);

    std::cout << "  num_partial_color_sets_very_dense = " << num_partial_color_sets_very_dense
              << " / " << num_total_partial_colors << " ("
              << (num_partial_color_sets_very_dense * 100.0) / num_total_partial_colors << "%)"
              << std::endl;
    std::cout << "  num_partial_color_sets_dense = " << num_partial_color_sets_dense << " / "
              << num_total_partial_colors << " ("
              << (num_partial_color_sets_dense * 100.0) / num_total_partial_colors << "%)"
              << std::endl;
    std::cout << "  num_partial_color_sets_sparse = " << num_partial_color_sets_sparse << " / "
              << num_total_partial_colors << " ("
              << (num_partial_color_sets_sparse * 100.0) / num_total_partial_colors << "%)"
              << std::endl;

    std::cout << "  partial colors: " << num_bits_colors / 8 << " bytes ("
              << (num_bits_colors * 100.0) / num_bits() << "%)\n";
    std::cout << "  meta colors: "
              << m_meta_color_sets.num_bytes() + m_meta_color_sets_offsets.num_bytes() << " bytes ("
              << ((m_meta_color_sets.num_bytes() + m_meta_color_sets_offsets.num_bytes()) * 8 *
                  100.0) /
                     num_bits()
              << "%)\n";
    std::cout << "  other: " << essentials::vec_bytes(m_partition_endpoints) << " bytes ("
              << ((essentials::vec_bytes(m_partition_endpoints) * 8) * 100.0) / num_bits()
              << "%)\n";
    std::cout << "  partition endpoints: ";
    for (auto p : m_partition_endpoints) std::cout << p.min_color << " ";
    std::cout << std::endl;
}

void differential::print_stats() const  //
{
    std::cout << "Color sets statistics:\n";
    std::cout << "  Number of partitions: " << num_partitions() << std::endl;

    uint64_t num_bits_representative_offsets = m_representative_offsets.num_bytes() * 8;
    uint64_t num_bits_color_sets_offsets = m_color_set_offsets.num_bytes() * 8;
    uint64_t num_bits_color_sets = m_color_sets.num_bytes() * 8;

    const uint64_t num_clusters = m_clusters.num_bits();
    uint64_t num_representatives = 0;
    uint64_t num_differential_color_sets = 0;
    uint64_t num_metadata = 0;

    uint64_t size_representatives = 0;
    uint64_t size_differentials = 0;

    uint64_t num_colors_tenth = num_colors() / 10;

    std::vector<uint64_t> distribution(11, 0);

    for (uint64_t representative_id = 0; representative_id < num_partitions();
         representative_id++)  //
    {
        uint64_t representative_begin = m_representative_offsets.access(representative_id);
        auto it = m_color_sets.get_iterator_at(representative_begin);
        uint64_t prev_position = it.position();

        uint64_t size = bits::util::read_delta(it);
        size_representatives += size;
        num_metadata += it.position() - prev_position;
        prev_position = it.position();

        for (uint64_t i = 0; i < size; i++) {
            bits::util::read_delta(it);
            num_representatives += it.position() - prev_position;
            prev_position = it.position();
        }
    }
    for (uint64_t color_id = 0; color_id < num_color_sets(); color_id++)  //
    {
        uint64_t color_set_begin = m_color_set_offsets.access(color_id);
        auto it = m_color_sets.get_iterator_at(color_set_begin);
        uint64_t prev_position = it.position();

        uint64_t size = bits::util::read_delta(it);
        size_differentials += size;
        num_metadata += it.position() - prev_position;
        prev_position = it.position();

        bits::util::read_delta(it);  // original color_set size
        num_metadata += it.position() - prev_position;
        prev_position = it.position();

        for (uint64_t i = 0; i < size; i++) {
            bits::util::read_delta(it);
            uint64_t delta_size = it.position() - prev_position;
            num_differential_color_sets += delta_size;

            prev_position = it.position();
        }
        uint64_t q = 0;
        if (num_colors_tenth != 0) {
            q = size / (num_colors_tenth) > 10 ? 10 : size / (num_colors_tenth);
        }

        distribution[q]++;
    }

    assert(num_bits() > 0);
    assert(num_bits_color_sets > 0);

    std::cout << "  representative offsets: " << num_bits_representative_offsets / 8 << " bytes ("
              << (num_bits_representative_offsets * 100.0) / num_bits() << "%)" << std::endl;
    std::cout << "  average representative set size: "
              << size_representatives * 1. / num_partitions() << " ints" << std::endl;
    std::cout << "  average differential set size: " << size_differentials * 1. / num_color_sets()
              << " ints" << std::endl;
    std::cout << "  differential color set offsets: " << num_bits_color_sets_offsets / 8
              << " bytes (" << (num_bits_color_sets_offsets * 100.0) / num_bits() << "%)"
              << std::endl;
    std::cout << "  clusters: " << num_clusters / 8 << " bytes ("
              << (num_clusters * 100.0) / num_bits() << "%)" << std::endl;
    std::cout << "  differential color sets: " << num_bits_color_sets / 8 << " bytes ("
              << (num_bits_color_sets * 100.0) / num_bits() << "%)" << std::endl;
    std::cout << "    representatives: " << num_representatives / 8 << " bytes ("
              << (num_representatives * 100.0) / num_bits_color_sets << "%)" << std::endl;
    std::cout << "    differential color sets: " << num_differential_color_sets / 8 << " bytes ("
              << (num_differential_color_sets * 100.0) / num_bits_color_sets << "%)" << std::endl;
    std::cout << "    metadata: " << num_metadata / 8 << " bytes ("
              << (num_metadata * 100.0) / num_bits_color_sets << "%)" << std::endl;
    std::cout << "  differential color sets size distribution:" << std::endl;
    for (uint64_t partition = 0; partition < distribution.size(); partition++) {
        std::cout << distribution[partition] << " ";
    }
    std::cout << std::endl;
}

void meta_differential::print_stats() const  //
{
    std::cout << "Color sets statistics:\n";
    std::cout << "  Number of partitions: " << num_partitions() << '\n';
    std::cout << "  Number of partition sets: " << num_partition_sets() << '\n';

    uint64_t num_bits_meta_color_sets =
        8 * (m_relative_colors_offsets.num_bytes() + m_partition_sets_offsets.num_bytes() +
             m_relative_colors.num_bytes() + m_partition_sets.num_bytes() +
             m_partition_sets_partitions.num_bytes() +
             m_partition_sets_partitions_rank1_index.num_bytes());

    uint64_t num_bits_partial_color_sets = 0;
    for (auto const& pcs : m_partial_color_sets) num_bits_partial_color_sets += pcs.num_bits();

    assert(num_bits() > 0);
    std::cout << "  partial color sets: " << num_bits_partial_color_sets / 8 << " bytes ("
              << (num_bits_partial_color_sets * 100.0) / num_bits() << "%)\n";
    std::cout << "  meta color sets: " << num_bits_meta_color_sets / 8 << " bytes ("
              << (num_bits_meta_color_sets * 100.0) / num_bits() << "%)\n";
    std::cout << "  other: " << essentials::vec_bytes(m_partition_endpoints) << " bytes ("
              << ((essentials::vec_bytes(m_partition_endpoints) * 8) * 100.0) / num_bits()
              << "%)\n";
}

}  // namespace fulgor
