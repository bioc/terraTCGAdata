.isSingleChar <- function(ch) {
    stopifnot(is.character(ch), !is.na(ch), identical(length(ch), 1L))
}

#' Obtain or set the Terra Workspace Project Dataset
#'
#' Terra allows access to about 71 open access TCGA datasets. A dataset
#' workspace can be set using the `terraTCGAworkspace` function with a
#' `projectName` input. Use the `selectTCGAworkspace` function to select
#' a TCGA data workspace from an interactive table.
#'
#' @details
#'   Note that GDC workspaces are not supported and are excluded
#'   from the search results. GDC workspaces use a Terra workflow to download
#'   TCGA data rather than providing Google Bucket storage locations for easy
#'   data retrieval. To reset the `terraTCGAworkspace`, use
#'   `terraTCGAworkspace(NULL)` and you will be prompted to select from a list
#'   of TCGA workspaces. You may also check the current active workspace by
#'   running `terraTCGAworkspace()` without any inputs.
#'
#' @aliases selectTCGAworkspace
#'
#' @param projectName character(1) A project code usually in the form of
#'   `TCGA_CODE_OpenAccess_V1-0_DATA`. See `selectTCGAworkspace` to
#'   interactively select from a table of project codes.
#'
#' @param verbose logical(1) Whether to provide more informative messages
#'   when an the "terraTCGAdata.workspace" option is set.
#'
#' @param ... further arguments passed down to lower level functions, not
#'   intended for the end user.
#'
#' @return A Terra TCGA Workspace name
#'
#' @md
#'
#' @examples
#' if (
#'     GCPtools::gcloud_exists() &&
#'     identical(AnVILBase::avplatform_namespace(), "AnVILGCP") &&
#'     nzchar(AnVILGCP::avworkspace_name())
#' ) {
#'   selectTCGAworkspace()
#'   terraTCGAworkspace()
#' }
#'
#' @export
terraTCGAworkspace <-
    function(projectName = getOption("terraTCGAdata.workspace", NULL))
{
    if (!is.null(projectName))
        opt <- .validateWorkspace(projectName = projectName)
    else
        opt <- selectTCGAworkspace(projectName = projectName)
    if (is.null(opt) || !nzchar(opt))
        warning("'terraTCGAdata.workspace' is blank; see '?terraTCGAworkspace'")
    opt
}

.validateWorkspace <- 
    function(projectName)
{
    .isSingleChar(projectName)
    tcga_choices <- .getWorkspaceTable()[["name"]]
    validPC <- projectName %in% tcga_choices
    if (!validPC)
        stop("'projectName' not in the 'selectTCGAworkspace()' table ")
    ws <- getOption("terraTCGAdata.workspace", NULL)
    if (!identical(ws, projectName) && !is.null(ws)) {
        warning(
            "'terraTCGAData.workspace' option set to ",
            sQuote(projectName, q = FALSE),
            " from ",
            sQuote(ws, q = FALSE),
            call. = FALSE
        )
        ws <- projectName
        options("terraTCGAdata.workspace" = projectName)
    }
    ws
}

.done_fun <- function(avs, row_selected) {
    unname(unlist(avs[row_selected, "name"]))
}

## from AnVIL:::.workspaces()
#' @importFrom AnVILGCP avworkspaces
.workspaces <- local({
    workspaces <- NULL
    function() {
        if (is.null(workspaces))
            workspaces <<- avworkspaces()
        workspaces
    }
})

#' Obtain the table of datasets from the Terra platform
#' 
#' The datasets include all TCGA datasets that do not come from the Genomic
#' Data Commons Data Repository because those data use a different data model.
#' 
#' @param project character(1) A prefix for the regex search across all public
#' projects on the terra platform (default: `"^TCGA"`). Usually, this does not
#' change.
#'
#' @param cancerCode character(1) Corresponds to the TCGA cancer code (e.g,
#'   "ACC" for AdrenoCortical Carcinoma) of interest. The default value of
#'   (`.*`) provides all available cancer datasets.
#'
#' @keywords internal
#' 
#' @return A `tibble` `data.frame` that match the project in put; by default,
#'   TCGA workspaces.
.getWorkspaceTable <- function(project = "^TCGA", cancerCode = ".*") {
    avs <- .workspaces()
    gdcind <- grep("GDCDR", avs[["name"]], invert = TRUE)
    avs <- avs[gdcind, ]
    project_code <- paste(project = project, cancerCode = cancerCode, sep = "_")
    ind <- grep(project_code, avs[["name"]])
    avs[ind, ]
}

#' @describeIn terraTCGAworkspace Function to interactively select from the
#'   available TCGA data workspaces in Terra. The `projectName` argument and
#'   'terraTCGAdata.workspace' option must be `NULL` to enable the interactive
#'   gadget.
#'
#' @export
selectTCGAworkspace <- function(
    projectName = getOption("terraTCGAdata.workspace", NULL), verbose = FALSE, ...
) {
    wst <- .getWorkspaceTable(...)
    if (interactive() && is.null(projectName)) {
        ws <- AnVIL::.gadget_run(
            "Terra TCGA Workspaces", wst, .done_fun
        )
        if (verbose)
            message("Setting 'terraTCGAdata.workspace' option to ", ws)
    } else {
        ws <- projectName    
    }
    options("terraTCGAdata.workspace" = ws)
    ws
}
