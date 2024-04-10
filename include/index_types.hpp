#pragma once

#include "index.hpp"

#include "builder.hpp"
#include "color_classes/color_classes.hpp"

namespace fulgor {
typedef index<hybrid> hybrid_colors_index_type;
typedef hybrid_colors_index_type index_type;  // in use
}  // namespace fulgor

#include "meta_builder.hpp"
#include "color_classes/meta.hpp"

namespace fulgor {
typedef index<meta<hybrid>> meta_hybrid_colors_index_type;
typedef meta_hybrid_colors_index_type meta_index_type;  // in use
}  // namespace fulgor

#include "differential_builder.hpp"
#include "color_classes/differential.hpp"

namespace fulgor {
typedef index<differential> differential_hybrid_colors_index_type;
typedef differential_hybrid_colors_index_type differential_index_type;  // in use
}