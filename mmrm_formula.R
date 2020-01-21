#' Mix model with repeated measurements (MMRM) model
#'
#' The MMRM table function summarizes MMRM test results by visit and groups. The
#' function produces adjusted lsmeans and standard error, as well as conducts
#' comparisons between groups' adjusted means, where the first level of the group
#' is the reference level.
#'
#' @details \url{https://pages.github.roche.com/NEST/docs/bookdown/NEST/statistics_clinical_trials/master/mixed-effect-model-repeat-measurement.html}
#'
#' Create rtable of MMRM test results
#'
#'
#'
#' @inheritParams a_mmrm
#' @param data a \code{data.frame} with all the variable that are used in
#'   \code{formula}. The \code{\link{arm}} and \code{\link{visit}} columns in data must be factors.
#' @param id_var a character describing the variable name used as subject IDs,
#'   "USUBJID" by default.
#' @param col_N vector that specifies the denominator (i.e. total number of
#'   subjects) in each \code{\link{arm}}.
#' @param mode algorithm for degree of freedom: "auto", "df.error" or
#'   "boot-satterthwaite".
#' @param conf.level confidence level. Must be number greater than 0 and less
#'   than 1.
#' @param weights argument from emmeans, "proportional" by default.
#' @param table_tree boolean to specify whether to output rtable as table_tree
#'
#' @return rtable object or table tree, depending on the the `table_tree` argument
#'
#' @export
#'
#' @seealso \code{\link{a_mmrm}}
#'
#' @examples
#' library(dplyr)
#' library(nlme)
#' library(emmeans)
#' library(random.cdisc.data)
#'
#' ADSL <- radsl(cached = TRUE)
#' ADQS <- radqs(cached = TRUE)
#' ADQS_f <- ADQS %>%
#'   filter(PARAMCD=="FKSI-FWB" & ABLFL != "Y" & ABLFL2 != "Y") %>%
#'   droplevels()
#'
#' t_mmrm(formula = AVAL ~ arm(ARM) + visit(AVISIT) + STRATA1 + BMRKR2 + BASE + ARM*AVISIT,
#'        data = ADQS_f,
#'        id_var = "USUBJID",
#'        col_N = table(ADSL$ARM),
#'        mode = "boot-satterthwaite",
#'        conf.level = 0.95,
#'        table_tree = FALSE,
#'        weights = "proportional",
#'        corStruct = "corSymm"
#' )
#'
t_mmrm <- function(formula = AVAL ~ arm(ARM) + visit(AVISIT) + ARM * VISIT,
                   data,
                   id_var = "USUBJID",
                   col_N = table(ADSL$ARM),
                   mode = c("df.error", "auto", "boot-satterthwaite"),
                   conf.level = 0.95,
                   weights = "proportional",
                   corStruct = NULL,
                   table_tree = TRUE) {

  mmrm_result <- a_mmrm(formula = formula,
                      data = data,
                      id_var = id_var,
                      mode = mode,
                      conf.level = conf.level,
                      weights = weights,
                      corStruct = corStruct)

  arm_var <- mmrm_result$`arm_var`
  visit_var <- mmrm_result$`visit_var`

  s_contrast_df <- split(
    mmrm_result$`contrast`,
    mmrm_result$`contrast`[visit_var]
  )

  s_estimate_df <- split(
    mmrm_result$`estimate`,
    mmrm_result$`estimate`[visit_var]
  )

  arm_lvl <- levels(data[[arm_var]])

  tbl_head <- rheader(rrowl("", arm_lvl))

  format_pval <- function(x, output) {
    if (x < 0.0001) {
      "<.0001"
    } else {
      paste(round(x, 4))
    }
  }

  mmrm_node_list <- Map(function(df1_i, df2_i, visit){

    tbl <- rtable(
      header = tbl_head,
      rrowl("n",
            tapply(df1_i$`n`, factor(df1_i[[arm_var]], levels = arm_lvl),
                   function(n_i) {
                     c(n_i)
                   }),
            format = "xx"
      ),
      rrowl("Adjusted Mean (SE)",
            lapply(split(df1_i,
                         factor(df1_i[[arm_var]], levels = arm_lvl),
                         drop = FALSE),
                   function(vector_i) {
                     if (is.null(vector_i)) {
                       NULL}
                     else {
                       c(vector_i$`emmean`, vector_i$`SE`)
                     }
                   }),
            format = sprintf_format("%.3f (%.3f)")
      ),
      rrowl(paste0(mmrm_result$`conf_level`*100, "% CI"),
            lapply(split(df1_i,
                         factor(df1_i[[arm_var]], levels = arm_lvl),
                         drop = FALSE),
                   function(vector_i) {
                     if (is.null(vector_i)) {
                       NULL}
                     else {
                       c(vector_i$`lower.CL`, vector_i$`upper.CL`)
                     }
                   }),
            format = "(xx.xxx, xx.xxx)"
      ),
      rrow(),
      rrowl(paste0("Difference in Adjusted Means (SE) (vs. ", mmrm_result$`ref_level`, ")"),
            lapply(split(df2_i,
                         factor(df2_i[[arm_var]], levels = arm_lvl),
                         drop = FALSE),
                   function(vector_i) {
                     if (is.null(vector_i)) {
                       NULL
                     } else {
                       c(vector_i$`estimate`, vector_i$`SE`)
                     }
                   }),
            format = sprintf_format("%.3f (%.3f)")
      ),
      rrowl(paste0(mmrm_result$`conf_level`*100, "% CI"),
            lapply(split(df2_i,
                         factor(df2_i[[arm_var]], levels = arm_lvl),
                         drop = FALSE),
                   function(vector_i) {
                     if (is.null(vector_i)) {
                       NULL
                     } else {
                       c(vector_i$`lower.CL`, vector_i$`upper.CL`)
                     }
                   }),
            format = "(xx.xxx, xx.xxx)"
      ),
      rrowl("Relative Reduction (%)",
            lapply(split(df2_i, df2_i[[arm_var]], drop = FALSE),
                   function(vector_i) {
                     if (is.null(vector_i)) {
                       NULL
                     } else {
                       c(vector_i$`relative_reduc`)
                     }
                   }),
            format = "xx.x%"
      ),
      rrow(),
      rrowl("p-value (MMRM)",
            lapply(split(df2_i, df2_i[[arm_var]], drop = FALSE),
                   function(vector_i) {
                     if (is.null(vector_i)) {
                       NULL
                     } else {
                       c(vector_i$`p.value`)
                     }
                   }),
            format = format_pval
      )
    )

    tbl <- header_add_N(tbl, col_N)

    node(
      name = visit,
      content = tbl,
      children = NULL
    )

  }, s_estimate_df, s_contrast_df, names(s_estimate_df))

  tree <- invisible_node(
    name = "root",
    children = mmrm_node_list,
    content = NULL
  )

  if (table_tree) {
    tree
  } else {
    to_rtable(tree)
  }
}

#' MMRM model, test, estimate
#'
#'
#' @param formula a gls formula, the arm and visit variable needs to be wrapped
#'   in \code{\link{arm}} and \code{\link{visit}}.
#' @param data a \code{data.frame} with all the variable that are used in
#'   \code{formula}. The arm column in data must be a factor.
#' @param id_var a character describing the variable name used as subject IDs,
#'   "USUBJID" by default.
#' @param mode algorithm for degree of freedom: "auto", "df.error" or
#'   "boot-satterthwaite".
#' @param conf.level confidence level. Must be number greater than 0 and less
#'   than 1.
#' @param weights argument from emmeans, "proportional" by default.
#' @param corStruct \code{NULL} by default or a string with the name of \link[nlme]{corClasses}.
#'
#' @return a dataframe with MMRM results
#'
#' @importFrom nlme gls corSymm corAR1 corCompSymm varIdent
#' @importFrom emmeans emmeans contrast
#'
#' @export
#'
#' @examples
#' library(random.cdisc.data)
#' library(dplyr)
#'
#' ADSL <- radsl(cached = TRUE)
#' ADQS <- radqs(cached = TRUE)
#' ADQS_f <- ADQS %>%
#'   filter(PARAMCD == "FKSI-FWB" & ABLFL != "Y" & ABLFL2 != "Y") %>%
#'   droplevels()
#'
#' mmrm_results <- a_mmrm(
#'   data = ADQS_f,
#'   formula = AVAL ~ arm(ARM) + visit(AVISIT) + STRATA1 + BMRKR2 + BASE + ARM * AVISIT,
#'   id_var = "USUBJID",
#'   mode = "boot-satterthwaite",
#'   conf.level = 0.95,
#'   weights = "proportional",
#'   corStruct = "corSymm"
#' )
#'
#' names(mmrm_results)
#'
#' mmrm_results["contrast"]
#' mmrm_results["estimate"]
#'
a_mmrm <- function(data,
                 formula = AVAL ~ arm(ARM) + visit(AVISIT) + ARM * VISIT,
                 id_var = "USUBJID",
                 mode = c("df.error", "auto", "boot-satterthwaite"),
                 conf.level = 0.95,
                 weights = "proportional",
                 corStruct = NULL
                 ) {

  mode <- match.arg(mode)

  stopifnot(
    is.data.frame(data),
    id_var %in% names(data),
    is.null(corStruct) || corStruct %in% c(
      "corAR1", "corARMA", "corCAR1", "corCompSymm", "corExp", "corGaus",
      "corLin", "corRatio", "corSpher", "corSymm"
    )
  )

  tm <- mmrm_items(formula = formula, data = data)

  form <- tm$formula

  # names are strings and symbols can be substituted into bquote and used with !!
  response_name <- tm$rsp_name
  arm_name <- tm$arm_name
  visit_name <- tm$visit_name

  arm_symbol <- sym(arm_name)

  arm_values <- data[[arm_name]]
  response_values <- data[[response_name]]

  stopifnot(nlevels(arm_values) >= 2)
  reference_level <- levels(arm_values)[1]

  stopifnot(
    is.numeric(response_values),
    is.factor(arm_values),
    is.factor(data[[visit_name]])
  )

  # Modeling step in SAS assumes complete case; however estimation step considers all data
  if (any(is.na(response_values))) {
    warning(
      "Some records have a missing endpoint, which will be excluded from MMRM modeling.",
      head(data[is.na(response_values), c(id_var, visit_name, response_name)]),
      call. = FALSE
    )
  }

  # remove all entries where response is NA, droplevels as well
  data_cc <- data %>%
    filter(!is.na(!!sym(response_name))) %>%
    droplevels()

  # check all arms will still be present after NA filtering
  stopifnot(nlevels(arm_values) == nlevels(data_cc[[arm_name]]))

  # each arm should have at least have 5 records
  if (!all(table(data_cc[[arm_name]]) > 5)) {
    stop(paste("Each group / arm should have at least 5 records with non-missing", response_name))
  }

  cor_formula <- as.formula(paste("~", "as.numeric(", visit_name, ") | ", id_var))

  correlation <- if (is.null(corStruct)) {
    NULL
  } else {
    do.call(corStruct, list(form = cor_formula))
  }

  fit_mmrm <- gls(
    model = form,
    correlation = correlation,
    weights = varIdent(form = as.formula(paste(" ~ 1 |", visit_name))),
    method = "REML",
    data = data_cc, # model fit on complete case
    na.action = na.exclude
  )

  # hacky for emmeans since emmeans evaluation and namespace is not good
  # https://github.com/rvlenth/emmeans/issues/13
  fit_mmrm$call$model <- substitute(form)
  emm <- emmeans(
    fit_mmrm,
    mode = mode,
    specs = as.formula(paste("~ ", arm_name, "|", visit_name)),
    data = data %>% select(-c(response_name)), # estimate on original data
    weights = weights
  )

  # Relative Reduction (in change from baseline) is calculated using model based
  # LS means as 100*(LS mean change from baseline in Control Pooled group –
  # LS mean change from baseline in Treatment Group)/LS mean change from
  # baseline in Control Pooled group.

  ### adjusted estimate for each arm
  estimate <- confint(emm, level = conf.level) %>%
    as.data.frame()

  # add n information

  data_n <- data %>%
    group_by_at(.vars = c(visit_name, arm_name)) %>%
    summarise(n = n())

  estimate <- estimate %>%
    left_join(., data_n, by = c(visit_name, arm_name))

  # get emmean for reference group to join into full dataframe so that relative reduction in
  # emmean (mean of response variable) can be computed with respect to reference level (e.g. ARM A)
  means_at_ref <- estimate %>%
    filter(!!arm_symbol == reference_level) %>%
    select(c(visit_name, "emmean")) %>%
    rename(ref = emmean)

  relative_reduc <- estimate %>%
    filter(!!arm_name != reference_level) %>%
    left_join(means_at_ref, by = c(visit_name)) %>%
    mutate(relative_reduc = (ref - emmean) / ref) %>%
    select(c(visit_name, arm_name, "relative_reduc"))

  sum_fit_diff <- summary(contrast(emm, method = "trt.vs.ctrl"), level = conf.level, infer = c(TRUE, TRUE), adjust = "none")

  # get the comparison group name from "contrast" column, e.g. "ARMB - ARMA" returns "ARMB", i.e. remove " - ARMA"
  contrast <- sum_fit_diff %>%
    mutate(
      col_by = gsub(paste0("\\s-\\s", reference_level), "", contrast) %>%
        factor(., levels = levels(data[[arm_name]]))
    ) %>%
    select(-contrast) %>%
    rename(!!arm_symbol := col_by) %>%
    left_join(., relative_reduc, by = c(visit_name, arm_name))

  list(
    contrast = contrast,
    estimate = estimate,
    arm_var = arm_name,
    visit_var = visit_name,
    ref_level = reference_level,
    conf_level = conf.level
  )

}

# the response, arm and visit variables are dealt with here
mmrm_items <- function(formula, data) {
  # Process the special variables in formula
  mt <- terms(formula, specials = c("arm", "visit"), data = data)

  if (!all(all.vars(attr(mt, "variables")) %in% names(data))) {
    stop("All formula variables must appear in 'data'")
  }
  # get indexes
  irsp <- attr(mt, "response")
  iarm <- attr(mt, "specials")$arm
  ivisit <- attr(mt, "specials")$visit

  stopifnot(
    !is.null(irsp),
    !is.null(ivisit),
    !is.null(iarm)
  )

  # process the special variables in formula
  arm_name <- all.vars(formula)[[iarm]]
  visit_name <- all.vars(formula)[[ivisit]]

  # update right hand side of formula (lhs ~ rhs), i.e. arm(ARM) -> ARM and interaction terms
  cleaned_formula <- paste(c(arm_name, visit_name, attr(mt, "term.labels")[-c(iarm - 1, ivisit - 1)]), collapse = "+") %>%
    paste(". ~", .) %>%
    update.formula(
      old = as.formula(mt),
      new = as.formula(.)
    )

  environment(cleaned_formula) <- new.env() # no scoping for formula elements needed

  list(
    formula = cleaned_formula,
    rsp_name = all.vars(formula)[irsp],
    arm_name = arm_name,
    visit_name = visit_name,
    model_frame = stats::model.frame(formula = cleaned_formula, data = data, na.action = stats::na.pass)
  )
}
