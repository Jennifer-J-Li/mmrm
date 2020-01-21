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
