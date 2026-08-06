#pragma once

#include "GenericNeighborSearch.h"
#include "alg/orders/RowMajorOrder.h"
#include "alg/orders/MortonOrder.h"

using RowMajorNeighborSearch = GenericNeighborSearch<RowMajorOrder>;
using MortonNeighborSearch = GenericNeighborSearch<MortonOrder>;
