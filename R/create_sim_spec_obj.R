create_simulatr_specifier_object <- function(param_grid, fixed_params, one_rep_times, methods = c("glmeiv_slow", "glmeiv_fast", "thresholding")) {
  methods <- sort(methods)

  ####################################
  # 1. Define data_generator function.
  ####################################
  data_generator_object <- simulatr::simulatr_function(f = generate_full_data,
                                                       arg_names = c("m_fam", "m_intercept", "m_perturbation", "g_fam", "g_intercept", "g_perturbation", "pi",
                                                                     "n", "B", "covariate_matrix", "m_covariate_coefs", "g_covariate_coefs", "m_offset", "g_offset",
                                                                     "run_unknown_theta_precomputation"),
                                                       packages = "glmeiv",
                                                       loop = FALSE,
                                                       one_rep_time = one_rep_times[["generate_data_function"]])

  ###########################
  # 2. Define the method list.
  ###########################
  method_list <- c(
    if ("glmeiv_fast" %in% methods) simulatr::simulatr_function(f = run_glmeiv_at_scale_simulatr,
                                                                arg_names = c("m_fam", "g_fam", "covariate_matrix", "m_offset", "g_offset", "alpha",
                                                                              "n_em_rep", "save_membership_probs_mult", "pi_guess_range",
                                                                              "m_perturbation_guess_range", "g_perturbation_guess_range"),
                                                                packages = "glmeiv",
                                                                loop = TRUE,
                                                                one_rep_time = one_rep_times[["glmeiv_fast"]]) else NULL,
    if ("glmeiv_slow" %in% methods) simulatr::simulatr_function(f = run_glmeiv_random_init_simulatr,
                                                                arg_names = c("m_fam", "g_fam", "covariate_matrix", "m_offset",
                                                                              "g_offset", "alpha", "n_em_rep", "save_membership_probs_mult",
                                                                              "pi_guess_range", "m_perturbation_guess_range", "g_perturbation_guess_range",
                                                                              "m_intercept_guess_range", "g_intercept_guess_range", "m_covariate_coefs_guess_range",
                                                                              "g_covariate_coefs_guess_range"),
                                                                packages = "glmeiv",
                                                                loop = TRUE,
                                                                one_rep_time = one_rep_times[["glmeiv_slow"]]) else NULL,
    if ("thresholding" %in% methods) simulatr::simulatr_function(f = run_thresholding_method_simulatr,
                                                                 arg_names = c("g_intercept", "g_perturbation",
                                                                               "g_fam", "m_fam", "pi", "covariate_matrix",
                                                                               "g_covariate_coefs", "m_offset", "g_offset", "alpha"),
                                                                 packages = "glmeiv", loop = TRUE, one_rep_time = one_rep_times[["thresholding"]]) else NULL
  )

}
