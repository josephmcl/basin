#pragma once 

#include "components.h"

#include "mkl.h"
#include <vector>

using namespace sbp_sat_2d;

auto compute_lambda_matrix(/*
        double     **lambdaA,
  const double     **F,
  const double     **MinvFT,
  const std::size_t n_batch,
  const components &sbp*/) -> void;

