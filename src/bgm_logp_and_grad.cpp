#include <RcppArmadillo.h>
#include "bgm_helper.h"
#include "bgm_logp_and_grad.h"
#include "common_helpers.h"
#include "explog_switch.h"



/**
 * Computes the log-pseudoposterior contribution for a single main-effect parameter (bgm model).
 *
 * For the specified variable, this function evaluates the log-pseudoposterior of either:
 *  - an ordinal threshold parameter (category-specific), or
 *  - a Blume–Capel main-effect parameter (linear or quadratic).
 *
 * The log-pseudoposterior combines:
 *  - Prior contribution: Beta prior on the logistic scale.
 *  - Sufficient statistic contribution: from category counts or Blume–Capel statistics.
 *  - Likelihood contribution: vectorized across all persons using the residual matrix.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - residual_matrix: Matrix of residual scores (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category: Category counts per variable (used for ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - baseline_category: Reference category for Blume–Capel centering.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - main_alpha, main_beta: Prior hyperparameters for the Beta prior.
 *  - variable: Index of the variable under consideration.
 *  - category: Category index (ordinal variables only).
 *  - parameter: Parameter index (Blume–Capel only: 0 = linear, 1 = quadratic).
 *
 * Returns:
 *  - The log-pseudoposterior value for the specified parameter.
 *
 * Notes:
 *  - Exactly one of `category` or `parameter` is relevant depending on variable type.
 *  - Uses a numerically stable denominator with exponential bounding.
 *  - This function is used within Metropolis and gradient-based updates.
 */
double log_pseudoposterior_main_effects_component (
    const arma::mat& main_effects,
    const arma::mat& residual_matrix,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const int variable,
    const int category,
    const int parameter
) {
  const int num_persons = residual_matrix.n_rows;
  double log_posterior = 0.0;

  auto log_beta_prior = [&](double main_effect_param) {
    return main_effect_param * main_alpha - std::log1p (MY_EXP (main_effect_param)) * (main_alpha + main_beta);
  };

  const int num_cats = num_categories(variable);

  if (is_ordinal_variable(variable)) {
    // Prior contribution + sufficient statistic
    const double value = main_effects(variable, category);
    log_posterior += value * counts_per_category(category + 1, variable);
    log_posterior += log_beta_prior (value);

    // Vectorized likelihood contribution
    // For each person, we compute the unnormalized log-likelihood denominator:
    //   denom = exp (-bound) + sum_c exp (main_effect_param_c + (c+1) * residual_score - bound)
    // Where:
    //   - residual_score is the summed interaction score excluding the variable itself
    //   - bound = num_cats * residual_score (for numerical stability)
    //   - main_effect_param_c is the main_effect parameter for category c (0-based)
    arma::vec residual_score = residual_matrix.col (variable);                     // rest scores for all persons
    arma::vec bound = num_cats * residual_score;                                  // numerical bound vector
    arma::vec denom = ARMA_MY_EXP (-bound);                                      // initialize with base term
    arma::vec main_effect_param = main_effects.row (variable).cols (0, num_cats - 1).t ();   // main_effect parameters

    for (int cat = 0; cat < num_cats; cat++) {
      arma::vec exponent = main_effect_param(cat) + (cat + 1) * residual_score - bound;       // exponent per person
      denom += ARMA_MY_EXP (exponent);                                           // accumulate exp terms
    }

    // We then compute the total log-likelihood contribution as:
    //   log_posterior -= bound + log (denom), summed over all persons
    log_posterior -= arma::accu (bound + ARMA_MY_LOG (denom));                    // total contribution
  } else {
    const double value = main_effects(variable, parameter);
    const double linear_main_effect = main_effects(variable, 0);
    const double quadratic_main_effect = main_effects(variable, 1);
    const int ref = baseline_category(variable);

    // Prior contribution + sufficient statistic
    log_posterior += value * blume_capel_stats(parameter, variable);
    log_posterior += log_beta_prior(value);

    // Vectorized likelihood contribution
    // For each person, we compute the unnormalized log-likelihood denominator:
    //   denom = sum_c exp (θ_lin * c + θ_quad * (c - ref)^2 + c * residual_score - bound)
    // Where:
    //   - θ_lin, θ_quad are linear and quadratic main_effects
    //   - ref is the reference category (used for centering)
    //   - bound = num_cats * residual_score (stabilizes exponentials)
    arma::vec residual_score = residual_matrix.col(variable);                     // rest scores for all persons
    arma::vec bound = num_cats * residual_score;                                  // numerical bound vector
    arma::vec denom(num_persons, arma::fill::zeros);                          // initialize denominator
    for (int cat = 0; cat <= num_cats; cat++) {
      int centered = cat - ref;                                               // centered category
      double quad_term = quadratic_main_effect * centered * centered;                    // precompute quadratic term
      double lin_term = linear_main_effect * cat;                                      // precompute linear term

      arma::vec exponent = lin_term + quad_term + cat * residual_score - bound;
      denom += ARMA_MY_EXP (exponent);                                           // accumulate over categories
    }

    // The final log-likelihood contribution is then:
    //   log_posterior -= bound + log (denom), summed over all persons
    log_posterior -= arma::accu (bound + ARMA_MY_LOG (denom));                    // total contribution
  }

  return log_posterior;
}



/**
 * Computes the log-pseudoposterior contribution for a single pairwise interaction (bgm model).
 *
 * The contribution consists of:
 *  - Sufficient statistic term: interaction × pairwise count.
 *  - Likelihood term: summed over all observations, using either
 *    * ordinal thresholds, or
 *    * Blume–Capel quadratic/linear main effects.
 *  - Prior term: Cauchy prior on the interaction coefficient (if active).
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters.
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - pairwise_stats: Sufficient statistics for pairwise counts.
 *  - var1, var2: Indices of the variable pair being updated.
 *
 * Returns:
 *  - The log-pseudoposterior value for the specified interaction parameter.
 *
 * Notes:
 *  - Bounds are applied for numerical stability in exponential terms.
 *  - The function assumes that `pairwise_effects` is symmetric.
 *  - Used within Metropolis and gradient-based updates of pairwise effects.
 */
double log_pseudoposterior_interactions_component (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double pairwise_scale,
    const arma::imat& pairwise_stats,
    const int var1,
    const int var2
) {
  const int num_observations = observations.n_rows;

  double log_pseudo_posterior = 2.0 * pairwise_effects(var1, var2) * pairwise_stats(var1, var2);

  for (int var : {var1, var2}) {
    int num_categories_var = num_categories (var);

    // Compute rest score: contribution from other variables
    arma::vec residual_scores = observations * pairwise_effects.col (var);
    arma::vec bounds = arma::max (residual_scores, arma::zeros<arma::vec> (num_observations)) * num_categories_var;
    arma::vec denominator = arma::zeros (num_observations);

    if (is_ordinal_variable (var)) {
      // Ordinal variable: denominator includes exp (-bounds)

      denominator += ARMA_MY_EXP (-bounds);
      for (int category = 0; category < num_categories_var; category++) {
        arma::vec exponent = main_effects (var, category) + (category + 1) * residual_scores - bounds;
        denominator += ARMA_MY_EXP(exponent);
      }

    } else {
      // Binary/categorical variable: quadratic + linear term
      const int ref_cat = baseline_category (var);
      for (int category = 0; category <= num_categories_var; category++) {
        int centered_cat = category - ref_cat;
        double lin_term = main_effects (var, 0) * category;
        double quad_term = main_effects (var, 1) * centered_cat * centered_cat;
        arma::vec exponent = lin_term + quad_term + category * residual_scores - bounds;
        denominator += ARMA_MY_EXP (exponent);
      }
    }

    // Subtract log partition function and bounds adjustment
    log_pseudo_posterior -= arma::accu (ARMA_MY_LOG (denominator));
    log_pseudo_posterior -= arma::accu (bounds);
  }

  // Add Cauchy prior terms for included pairwise effects
  if (inclusion_indicator (var1, var2) == 1) {
    log_pseudo_posterior += R::dcauchy (pairwise_effects (var1, var2), 0.0, pairwise_scale, true);
  }

  return log_pseudo_posterior;
}



/**
 * Computes the full log-pseudoposterior for the bgm model.
 *
 * The log-pseudoposterior combines:
 *  - Main-effect contributions:
 *    * Ordinal variables: one parameter per category with Beta prior.
 *    * Blume–Capel variables: linear and quadratic parameters with Beta priors.
 *  - Pairwise-effect contributions:
 *    * Included interactions (per inclusion_indicator) with Cauchy prior.
 *  - Likelihood contributions:
 *    * Vectorized over all persons, with numerically stabilized denominators.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - pairwise_effects: Symmetric matrix of pairwise interaction strengths.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category: Category counts per variable (for ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - main_alpha, main_beta: Hyperparameters for the Beta priors.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - pairwise_stats: Pairwise sufficient statistics.
 *  - residual_matrix: Matrix of residual scores (persons × variables).
 *
 * Returns:
 *  - The scalar log-pseudoposterior value for the full model.
 *
 * Notes:
 *  - Exponential terms are bounded with nonnegative `bound` values for stability.
 *  - Pairwise effects are included only when marked in `inclusion_indicator`.
 *  - This is the top-level function combining both main and interaction components.
 */
double log_pseudoposterior (
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const arma::imat& pairwise_stats,
    const arma::mat& residual_matrix
) {

  const int num_variables = observations.n_cols;
  const int num_persons = observations.n_rows;

  double log_pseudoposterior = 0.0;

  // Calculate the contribution from the data and the prior
  auto log_beta_prior = [&](double main_effect_param) {
    return main_effect_param * main_alpha - std::log1p (MY_EXP (main_effect_param)) * (main_alpha + main_beta);
  };

  for (int variable = 0; variable < num_variables; variable++) {
    if (is_ordinal_variable(variable)) {
      const int num_cats = num_categories(variable);
      for (int cat = 0; cat < num_cats; cat++) {
        double value = main_effects(variable, cat);
        log_pseudoposterior += counts_per_category(cat + 1, variable) * value;
        log_pseudoposterior += log_beta_prior(value);
      }
    } else {
      double value = main_effects(variable, 0);
      log_pseudoposterior += log_beta_prior(value);
      log_pseudoposterior += blume_capel_stats(0, variable) * value;

      value = main_effects(variable, 1);
      log_pseudoposterior += log_beta_prior(value);
      log_pseudoposterior += blume_capel_stats(1, variable) * value;
    }
  }
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      if (inclusion_indicator(var1, var2) == 0) continue;

      double value = pairwise_effects(var1, var2);
      log_pseudoposterior += 2.0 * pairwise_stats(var1, var2) * value;
      log_pseudoposterior += R::dcauchy(value, 0.0, pairwise_scale, true); // Cauchy prior
    }
  }

  // Calculate the log denominators
  for (int variable = 0; variable < num_variables; variable++) {
    const int num_cats = num_categories(variable);
    arma::vec residual_score = residual_matrix.col (variable);                    // rest scores for all persons
    arma::vec bound = num_cats * residual_score;                                  // numerical bound vector
    bound = arma::clamp(bound, 0.0, arma::datum::inf);                        // only positive bounds to prevent overflow

    arma::vec denom;
    if (is_ordinal_variable(variable)) {
      denom = ARMA_MY_EXP (-bound);                                     // initialize with base term
      arma::vec main_effect_param = main_effects.row (variable).cols (0, num_cats - 1).t ();   // main_effect parameters for variable
      for (int cat = 0; cat < num_cats; cat++) {
        arma::vec exponent = main_effect_param(cat) + (cat + 1) * residual_score - bound; // exponent per person
        denom += ARMA_MY_EXP (exponent);                                          // accumulate exp terms
      }
    } else {
      const double lin_effect = main_effects(variable, 0);
      const double quad_effect = main_effects(variable, 1);
      const int ref = baseline_category(variable);

      denom.zeros(num_persons);
      for (int cat = 0; cat <= num_cats; cat++) {
        int centered = cat - ref;                                               // centered category
        double quad = quad_effect * centered * centered;                    // precompute quadratic term
        double lin = lin_effect * cat;                                      // precompute linear term
        arma::vec exponent = lin + quad + cat * residual_score - bound;
        denom += ARMA_MY_EXP (exponent);                                           // accumulate over categories
      }
    }

    log_pseudoposterior -= arma::accu (bound + ARMA_MY_LOG (denom));                    // total contribution
  }

  return log_pseudoposterior;
}



/**
 * Computes the gradient of the log-pseudoposterior for main and active pairwise parameters.
 *
 * Gradient components:
 *  - Observed sufficient statistics (from counts_per_category, blume_capel_stats, pairwise_stats).
 *  - Minus expected sufficient statistics (computed via probabilities over categories).
 *  - Plus gradient contributions from priors:
 *    * Beta priors on main effects.
 *    * Cauchy priors on active pairwise effects.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - pairwise_effects: Symmetric matrix of pairwise interaction strengths.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category: Category counts per variable (for ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - main_alpha, main_beta: Hyperparameters for Beta priors.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - pairwise_stats: Sufficient statistics for pairwise effects.
 *  - residual_matrix: Matrix of residual scores (persons × variables).
 *
 * Returns:
 *  - A vector containing the gradient of the log-pseudoposterior with respect to
 *    all main and active pairwise parameters, in the same order as
 *    `vectorize_model_parameters_bgm()`.
 */
arma::vec gradient_log_pseudoposterior(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const arma::imat& pairwise_stats,
    const arma::mat& residual_matrix
) {
  const int num_variables = observations.n_cols;
  const int num_persons   = observations.n_rows;
  const int num_main      = count_num_main_effects(num_categories, is_ordinal_variable);
  arma::imat index_matrix(num_variables, num_variables, arma::fill::zeros);

  // Count active pairwise effects + Index map for pairwise parameters
  int num_active = 0;
  for (int i = 0; i < num_variables - 1; i++) {
    for (int j = i + 1; j < num_variables; j++) {
      if (inclusion_indicator(i, j) == 1){
        index_matrix(i, j) = num_main + num_active++;
        index_matrix(j, i) = index_matrix(i, j);
      }
    }
  }

  // Allocate gradient vector (main + active pairwise only)
  arma::vec gradient(num_main + num_active, arma::fill::zeros);

  // ---- STEP 1: Observed statistics ----
  int offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    if (is_ordinal_variable(variable)) {
      const int num_cats = num_categories(variable);
      for (int cat = 0; cat < num_cats; cat++) {
        gradient(offset + cat) = counts_per_category(cat + 1, variable);
      }
      offset += num_cats;
    } else {
      gradient(offset) = blume_capel_stats(0, variable);
      gradient(offset + 1) = blume_capel_stats(1, variable);
      offset += 2;
    }
  }
  for (int i = 0; i < num_variables - 1; i++) {
    for (int j = i + 1; j < num_variables; j++) {
      if (inclusion_indicator(i, j) == 0) continue;
      int location = index_matrix(i, j);
      gradient(location) = 2.0 * pairwise_stats(i, j);
    }
  }

  // ---- STEP 2: Expected statistics ----
  offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    const int num_cats = num_categories(variable);
    arma::vec residual_score = residual_matrix.col(variable);
    arma::vec bound = num_cats * residual_score;
    bound = arma::clamp(bound, 0.0, arma::datum::inf);

    if (is_ordinal_variable(variable)) {
      arma::vec main_param = main_effects.row(variable).cols(0, num_cats - 1).t();
      bound += main_param.max();

      arma::mat exponents(num_persons, num_cats);
      for (int cat = 0; cat < num_cats; cat++) {
        exponents.col(cat) = main_param(cat) + (cat + 1) * residual_score - bound;
      }

      arma::mat probs = ARMA_MY_EXP(exponents);
      arma::vec denom = arma::sum(probs, 1) + ARMA_MY_EXP(-bound);
      probs.each_col() /= denom;

      // main effects
      for (int cat = 0; cat < num_cats; cat++) {
        gradient(offset + cat) -= arma::accu(probs.col(cat));
      }

      // pairwise effects
      for (int j = 0; j < num_variables; j++) {
        if (inclusion_indicator(variable, j) == 0 || variable == j) continue;
        arma::vec expected_value = arma::zeros(num_persons);
        for (int cat = 0; cat < num_cats; cat++) {
          expected_value += (cat + 1) * probs.col(cat) % observations.col(j);
        }
        int location = (variable < j) ? index_matrix(variable, j) : index_matrix(j, variable);
        gradient(location) -= arma::accu(expected_value);
      }
      offset += num_cats;
    } else {
      const int ref = baseline_category(variable);
      const double lin_eff = main_effects(variable, 0);
      const double quad_eff = main_effects(variable, 1);

      arma::mat exponents(num_persons, num_cats + 1);
      for (int cat = 0; cat <= num_cats; cat++) {
        int score = cat;
        int centered = score - ref;
        double lin  = lin_eff * score;
        double quad = quad_eff * centered * centered;
        exponents.col(cat) = lin + quad + score * residual_score - bound;
      }
      arma::mat probs = ARMA_MY_EXP(exponents);
      arma::vec denom = arma::sum(probs, 1);
      probs.each_col() /= denom;

      arma::ivec lin_score  = arma::regspace<arma::ivec>(0, num_cats);
      arma::ivec quad_score = arma::square(lin_score - ref);

      // main effects
      gradient(offset) -= arma::accu(probs * lin_score);
      gradient(offset + 1) -= arma::accu(probs * quad_score);

      // pairwise effects
      for (int j = 0; j < num_variables; j++) {
        if (inclusion_indicator(variable, j) == 0 || variable == j) continue;
        arma::vec expected_value = arma::zeros(num_persons);
        for (int cat = 0; cat < num_cats; cat++) {
          expected_value += (cat + 1) * probs.col(cat + 1) % observations.col(j);
        }
        int location = (variable < j) ? index_matrix(variable, j) : index_matrix(j, variable);
        gradient(location) -= arma::accu(expected_value);
      }
      offset += 2;
    }
  }

  // ---- STEP 3: Priors ----
  offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    if (is_ordinal_variable(variable)) {
      const int num_cats = num_categories(variable);
      for (int cat = 0; cat < num_cats; cat++) {
        const double p = 1.0 / (1.0 + MY_EXP(-main_effects(variable, cat)));
        gradient(offset + cat) += main_alpha - (main_alpha + main_beta) * p;
      }
      offset += num_cats;
    } else {
      for (int k = 0; k < 2; k++) {
        const double param = main_effects(variable, k);
        const double p = 1.0 / (1.0 + MY_EXP(-param));
        gradient(offset + k) += main_alpha - (main_alpha + main_beta) * p;
      }
      offset += 2;
    }
  }
  for (int i = 0; i < num_variables - 1; i++) {
    for (int j = i + 1; j < num_variables; j++) {
      if (inclusion_indicator(i, j) == 0) continue;
      int location = index_matrix(i, j);
      const double effect = pairwise_effects(i, j);
      gradient(location) -= 2.0 * effect / (effect * effect + pairwise_scale * pairwise_scale);
    }
  }

  return gradient;
}



/**
 * Computes the log-likelihood ratio for updating a single variable’s parameter (bgm model).
 *
 * The ratio compares the likelihood of the current parameter state versus a proposed state,
 * given the observed data, main effects, and residual contributions from other variables.
 *
 * Calculation:
 *  - Removes the current interaction contribution from the residual scores.
 *  - Recomputes denominators of the softmax likelihood under both current and proposed states.
 *  - Returns the accumulated log difference across all persons.
 *
 * Inputs:
 *  - variable: Index of the variable being updated.
 *  - interacting_score: Integer vector of interaction scores with other variables.
 *  - proposed_state: Candidate value for the parameter being updated.
 *  - current_state: Current value of the parameter.
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - num_categories: Number of categories per variable.
 *  - residual_matrix: Matrix of residual scores (persons × variables).
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *
 * Returns:
 *  - The log-likelihood ratio (current − proposed) for the specified variable.
 *
 * Notes:
 *  - For ordinal variables, the likelihood includes thresholds and category-specific scores.
 *  - For Blume–Capel variables, linear and quadratic terms are applied with centered categories.
 *  - Bounds are used to stabilize exponentials in the softmax denominator.
 */
double compute_log_likelihood_ratio_for_variable (
    int variable,
    const arma::ivec& interacting_score,
    double proposed_state,
    double current_state,
    const arma::mat& main_effects,
    const arma::ivec& num_categories,
    const arma::mat& residual_matrix,
    const arma::imat& observations,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category
) {
  // Convert interaction score vector to double precision
  arma::vec interaction = arma::conv_to<arma::vec>::from (interacting_score);

  const int num_persons = residual_matrix.n_rows;
  const int num_categories_var = num_categories (variable);

  // Compute adjusted linear predictors without the current interaction
  arma::vec residual_scores = residual_matrix.col (variable) - interaction * current_state;

  // Stability bound for softmax (scaled by number of categories)
  arma::vec bounds = arma::max (residual_scores, arma::zeros<arma::vec> (num_persons)) * num_categories_var;

  arma::vec denom_current = arma::zeros (num_persons);
  arma::vec denom_proposed = arma::zeros (num_persons);

  if (is_ordinal_variable (variable)) {
    denom_current += ARMA_MY_EXP(-bounds);
    denom_proposed += ARMA_MY_EXP(-bounds);

    for (int category = 0; category < num_categories_var; category++) {
      const double main = main_effects(variable, category);
      const int score = category + 1;

      for (int person = 0; person < num_persons; person++) {
        const double base = main + score * residual_scores[person] - bounds[person];

        const double exp_current = MY_EXP(base + score * interaction[person] * current_state);
        const double exp_proposed = MY_EXP(base + score * interaction[person] * proposed_state);

        denom_current[person] += exp_current;
        denom_proposed[person] += exp_proposed;
      }
    }
  } else {
    // Binary or categorical variable: linear + quadratic score
    const int ref_cat = baseline_category (variable);

    for (int category = 0; category <= num_categories_var; category++) {
      int centered = category - ref_cat;
      double lin_term = main_effects (variable, 0) * category;
      double quad_term = main_effects (variable, 1) * centered * centered;
      arma::vec exponent = lin_term + quad_term + category * residual_scores - bounds;

      denom_current += ARMA_MY_EXP (exponent + category * interaction * current_state);
      denom_proposed += ARMA_MY_EXP (exponent + category * interaction * proposed_state);
    }
  }

  // Accumulated log-likelihood difference across persons
  return arma::accu (ARMA_MY_LOG (denom_current) - ARMA_MY_LOG (denom_proposed));
}



/**
 * Computes the log-pseudolikelihood ratio for updating a single pairwise interaction (bgm model).
 *
 * The ratio compares the pseudo-likelihood under a proposed value versus the current value
 * of an interaction parameter between two variables.
 *
 * Calculation:
 *  1. Direct contribution from the interaction term:
 *       Δβ × ∑(score_var1 × score_var2) over all persons,
 *     where Δβ = proposed_state − current_state.
 *  2. Change in pseudo-likelihood for variable1, accounting for its updated interaction with variable2.
 *  3. Symmetric change in pseudo-likelihood for variable2, accounting for its updated interaction with variable1.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of pairwise interaction parameters.
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - num_persons: Number of persons (rows of observations).
 *  - variable1, variable2: Indices of the variable pair being updated.
 *  - proposed_state: Candidate value for the interaction parameter.
 *  - current_state: Current value of the interaction parameter.
 *  - residual_matrix: Matrix of residual scores (persons × variables).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - pairwise_stats: Sufficient statistics for pairwise counts.
 *
 * Returns:
 *  - The log-pseudolikelihood ratio (current − proposed) for the specified interaction.
 *
 * Notes:
 *  - Calls `compute_log_likelihood_ratio_for_variable()` for both variables involved.
 *  - The factor of 2.0 in the direct term reflects symmetry of the pairwise sufficient statistic.
 *  - Used in Metropolis updates for pairwise interaction parameters.
 */
double log_pseudolikelihood_ratio_interaction (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const int num_persons,
    const int variable1,
    const int variable2,
    const double proposed_state,
    const double current_state,
    const arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& pairwise_stats
) {
  double log_ratio = 0.0;
  const double delta = proposed_state - current_state;

  // Extract score vectors for both variables across all persons
  arma::ivec score1 = observations.col(variable1);
  arma::ivec score2 = observations.col(variable2);

  // (1) Direct interaction contribution to the linear predictor:
  //     Δβ × ∑(score1_i × score2_i) for all persons i
  log_ratio += 2.0 * pairwise_stats(variable1, variable2) * delta;

  // (2) Change in pseudo-likelihood for variable1 due to the update in its interaction with variable2
  log_ratio += compute_log_likelihood_ratio_for_variable (
    variable1, score2, proposed_state, current_state, main_effects,
    num_categories, residual_matrix, observations, is_ordinal_variable,
    baseline_category
  );

  // (3) Symmetric change for variable2 due to its interaction with variable1
  log_ratio += compute_log_likelihood_ratio_for_variable (
    variable2, score1, proposed_state, current_state, main_effects,
    num_categories, residual_matrix, observations, is_ordinal_variable,
    baseline_category
  );

  return log_ratio;
}
