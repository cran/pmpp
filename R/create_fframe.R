#' Add empty rows with time stamps to each cress-sectional unit in the panel
#'
#' @author Michal Oleszak
#'
#' @param indata \code{data.frame} with a panel structure
#' @param timestamps \code{vector} of time IDs for the added time periods
#' @param panel_ind vector of length 2 indicating names of variables indexing
#'                   units and time periods respectively
#' @param overwrite logical; if TRUE, existing rows in the data are overwritten
#'                  with empty rows if their time ID is in timestamps
#'
#' @description Creates a forecast frame as required by the \code{predict.pmpp()} method.
#'              To each cross-sectional unit in the data, a specified
#'              number of rows are added that contain only this unit's ID
#'              and the selected time ID.
#'
#' @return A \code{data.frame} with empty rows added.
#'
#' @importFrom dplyr full_join arrange_
#' @export
#'
#' @examples
#' data(EmplUK, package = "plm")
#' EmplUK <- dplyr::filter(EmplUK, year %in% c(1978, 1979, 1980, 1981, 1982))
#' my_fframe <- create_fframe(EmplUK, 1983)
create_fframe <- function(indata, timestamps, panel_ind = colnames(indata[, 1:2]), overwrite = FALSE) {
  unit_id <- panel_ind[1]
  time_id <- panel_ind[2]
  rows_to_add <- expand.grid(unique(indata[[unit_id]]), timestamps)
  colnames(rows_to_add) <- panel_ind
  if (!overwrite) {
    already_in_indata <- which(paste0(rows_to_add[, 1], rows_to_add[, 2]) %in%
      paste0(indata[[unit_id]], indata[[time_id]]))
    if (length(already_in_indata) > 0) {
      rows_to_add <- rows_to_add[-already_in_indata, ]
    }
  }
  fframe <- indata %>%
    full_join(rows_to_add, by = c(unit_id, time_id)) %>%
    arrange_(unit_id, time_id)
  return(fframe)
}
