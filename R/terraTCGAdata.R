.DEFAULT_TABLENAME <- "sample"
.DEFAULT_NAMESPACE <- "broad-firecloud-tcga"

#' Obtain a reference table for assay data in the Terra data model
#'
#' The column names in the output can be used in the `getAssayData` function.
#'
#' @inheritParams getAssayData
#'
#' @md
#'
#' @return A tibble of pointers to resources within the Terra data model
#'
#' @examples
#' if (
#'     GCPtools::gcloud_exists() &&
#'     identical(AnVILBase::avplatform_namespace(), "AnVILGCP") &&
#'     nzchar(AnVILGCP::avworkspace_name())
#' )
#'   getAssayTable(workspace = "TCGA_COAD_OpenAccess_V1-0_DATA")
#'
#' @export
getAssayTable <-
    function(
        tablename = .DEFAULT_TABLENAME, metacols = .PARTICIPANT_METADATA_COLS,
        workspace = terraTCGAworkspace(), namespace = .DEFAULT_NAMESPACE
    )
{
    .isGDC(workspace)
    samples <- avtable(
        table = tablename, namespace = namespace, name = workspace
    )
    metadata <- samples[, names(samples) %in% metacols]
    samples <- samples[, !names(samples) %in% metacols]
    samples[, grep("clin", names(samples), invert = TRUE)]
}

#' Obtain assay datasets from Terra
#'
#' @inheritParams getClinical
#'
#' @param assayName character() The name of the assay dataset column from
#'     `getAssayTable` to import into the current workspace.
#'
#' @param sampleCode character(1) The sample code used to filtering samples
#'     e.g., "01" for Primary Solid Tumors, see
#'     `data("sampleTypes", package = "TCGAutils")` for reference
#'
#' @param sampleIdx numeric() index or TRUE. Specify an index for subsetting the
#'     assay data. This argument is mainly used for example and vignette
#'     purposes. To use all the data, use the default value (default: `TRUE`)
#'
#' @return Either a matrix or RaggedExperiment depending on the assay selected
#'
#' @seealso [getAssayTable()]
#'
#' @md
#'
#' @examples
#' if (
#'     GCPtools::gcloud_exists() &&
#'     identical(AnVILBase::avplatform_namespace(), "AnVILGCP") &&
#'     nzchar(AnVILGCP::avworkspace_name())
#' )
#'   getAssayData(
#'       assayName = "protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
#'       sampleCode = c("01", "10"),
#'       workspace = "TCGA_ACC_OpenAccess_V1-0_DATA"
#'   )
#'
#' @export
getAssayData <-
    function(assayName, sampleCode = "01", tablename = .DEFAULT_TABLENAME,
        workspace = terraTCGAworkspace(), namespace = .DEFAULT_NAMESPACE,
        metacols = .PARTICIPANT_METADATA_COLS, sampleIdx = TRUE
    )
{
    if (missing(assayName))
        stop("Select an assay name from 'getAssayTable'")

    assayTable <- getAssayTable(
        tablename = tablename, metacols = metacols,
        workspace = workspace, namespace = namespace
    )
    assaylinks <- stats::na.omit(assayTable[, assayName])
    assayfiles <- unlist(unique(assaylinks))
    if (length(assayfiles[sampleIdx]) > 800)
        assayfiles <- split(
            assayfiles,
            cut(seq_along(assayfiles), 3,
            labels = utils::head(letters, 3))
        )

    bfc <- BiocFileCache::BiocFileCache()
    rpath <- BiocFileCache::bfcquery(bfc, assayName, exact = TRUE)[["rpath"]]
    if (!length(rpath)) {
        rpath <- BiocFileCache::bfcnew(bfc, rname = assayName, rtype = "local")
        dir.create(rpath)
        if (is.list(assayfiles))
            lapply(assayfiles, avcopy, destination = rpath)
        else
            avcopy(assayfiles[sampleIdx], rpath)
    }

    assaydl <- list.files(rpath, pattern = "\\.txt", full.names = TRUE)
    dfiles <- assaydl[sampleIdx]
    suffix <-
        if (grepl("^cna|snp", assayName)) { ".seg.txt" } else { ".data.txt" }
    tcgaids <- gsub(suffix, "", basename(dfiles), fixed = TRUE)
    if (!missing(sampleCode) && .is_character(sampleCode) && !is.null(sampleCode))
        tcgaids <- tcgaids[
            TCGAutils::TCGAsampleSelect(tcgaids, sampleCodes = sampleCode)
        ]
    assayselect <- file.path(rpath, paste0(tcgaids, suffix))
    datarows <- lapply(
        assayselect,
        readr::read_tsv,
        show_col_types = FALSE
    )
    if (grepl("seg", suffix)) {
        .mergeRangedSamples(datarows)
    } else {
        cnames <- unlist(datarows[[1]][1, ], use.names = FALSE)
        datarows <- lapply(datarows, function(x) x[-1, ])
        rowlist <- lapply(datarows, `[[`, 1L)
        dataonly <- lapply(datarows, function(x) {
            suppressMessages(readr::type_convert(x[2L]))
        })
        dups <- rowlist[!duplicated(rowlist)]
        logilist <- vector("list", length(dups))
        bindlist <- vector("list", length(dups))
        for (i in seq_along(dups)) {
            logilist[[i]] <-
                vapply(
                    rowlist,
                    function(rl) identical(rl, dups[[i]]),
                    logical(1L)
                )
            bindlist[[i]] <-
                cbind(rownames = dups[[i]],
                    dplyr::bind_cols(dataonly[logilist[[i]]]))
        }
        df <- Reduce(function(...) {
            merge(..., all = TRUE, by = "rownames")
        }, bindlist)
        rownames(df) <- df[["rownames"]]
        data.matrix(df[, -1])
    }
}

#' Import Terra TCGA data as a list
#'
#' @inheritParams getAssayData
#'
#' @param assayNames character() A vector of assays selected from the colnames
#'     of `getAssayTable`.
#'
#' @param verbose logical(1L) Whether to output additional details of the
#'   data facilitation.
#'
#' @return A `list` of assay datasets
#'
#' @md
#'
#' @examples
#' if (
#'     GCPtools::gcloud_exists() &&
#'     identical(AnVILBase::avplatform_namespace(), "AnVILGCP") &&
#'     nzchar(AnVILGCP::avworkspace_name())
#' )
#'   getTCGAdatalist(
#'       assayNames = c("protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
#'       "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg"),
#'       sampleCode = c("01", "10"),
#'       workspace = "TCGA_COAD_OpenAccess_V1-0_DATA"
#'   )
#'
#' @export
getTCGAdatalist <-
    function(
        assayNames, sampleCode, workspace = terraTCGAworkspace(),
        namespace = .DEFAULT_NAMESPACE, tablename = .DEFAULT_TABLENAME,
        sampleIdx = TRUE, verbose = TRUE 
    )
{
    nwspace <- paste0(namespace, "/", workspace)
    if (verbose)
        message("Using namespace/workspace: ", nwspace)
    dlist <- structure(vector("list", length(assayNames)), .Names = assayNames)
    for (assay in assayNames) {
        dlist[[assay]] <-
            getAssayData(
                assay, sampleCode = sampleCode, tablename = tablename,
                workspace = workspace, namespace = namespace,
                sampleIdx = sampleIdx
            )
    }
    dlist
}

#' Obtain a MultiAssayExperiment from the Terra workspace
#'
#' Workspaces on Terra come pre-loaded with TCGA Data. The examples in the
#' documentation correspond to the TCGA_COAD_OpenAccess_V1 workspace that
#' can be found on \url{app.terra.bio}.
#'
#' @inheritParams getClinical
#' @inheritParams getAssayData
#'
#' @param clinicalName character(1) The column name taken from
#'     `getClinicalTable()` and downloaded to be included as the `colData`.
#'
#' @param assays character() A character vector of assay names taken from
#'     `getAssayTable()`
#'
#' @param sampleCode character() A character vector of sample codes from
#'     `sampleTypesTable()`. By default, (NULL) all samples are downloaded and
#'     kept in the data.
#'
#' @param split logical(1L) Whether or not to split the `MultiAssayExperiment`
#'     by sample types using `splitAssays` helper function (default FALSE).
#'
#' @return A `MultiAssayExperiment` object with n number of assays corresponding
#'     to the `assays` argument.
#'
#' @md
#'
#' @examples
#' if (
#'     GCPtools::gcloud_exists() &&
#'     identical(AnVILBase::avplatform_namespace(), "AnVILGCP")
#' )
#'   terraTCGAdata(
#'       clinicalName = "clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin",
#'       assays = c("protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
#'       "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data"),
#'       workspace = "TCGA_COAD_OpenAccess_V1-0_DATA",
#'       sampleCode = NULL,
#'       sampleIdx = 1:4,
#'       split = FALSE
#'   )
#'
#' @export
terraTCGAdata <-
    function(
        clinicalName, assays, participants = TRUE,
        sampleCode = NULL, split = FALSE,
        workspace = terraTCGAworkspace(), namespace = .DEFAULT_NAMESPACE,
        tablename = .DEFAULT_TABLENAME, verbose = TRUE, sampleIdx = TRUE
    )
{
    nwspace <- paste0(namespace, "/", workspace)
    if (verbose)
        message("Using namespace/workspace: ", nwspace)
    datalist <- getTCGAdatalist(
        assayNames = assays, sampleCode = sampleCode, workspace = workspace,
        namespace = namespace, tablename = tablename, sampleIdx = sampleIdx,
        verbose = verbose 
    )
    explist <- methods::as(datalist, "ExperimentList")
    partIds <- TCGAutils::TCGAbarcode(
        unique(unlist(colnames(explist), use.names = FALSE))
    )
    coldata <- getClinical(
        columnName = clinicalName, workspace = workspace,
        namespace = namespace, participantIds = partIds, verbose = verbose 
    )
    coldata <- .transform_clinical_to_coldata(
        clinical_data = coldata,
        columnName = clinicalName,
        participants = participants,
        workspace = workspace,
        namespace = namespace
    )
    samap <- TCGAutils::generateMap(datalist, coldata, TCGAutils::TCGAbarcode)
    mae <- MultiAssayExperiment(
        experiments = datalist, colData = coldata, sampleMap = samap
    )
    if (split && (length(sampleCode) || !is.null(sampleCode))) {
        TCGAutils::TCGAsplitAssays(
            mae, sampleCodes = sampleCode, exclusive = TRUE
        )
    } else {
        mae
    }
}

#' @importFrom AnVILGCP avtable
.transform_clinical_to_coldata <-
    function(clinical_data, columnName, participants, workspace, namespace)
{
    coldata <- as.data.frame(clinical_data)
    if (!is.null(coldata[["patient.bcr_patient_barcode"]])) {
        rownames(coldata) <- coldata[["participant_id"]] <-
            toupper(coldata[["patient.bcr_patient_barcode"]])
        if (participants) {
            parts <- avtable(
                table = "participant", namespace = namespace, name = workspace
            )
            parts <- methods::as(parts, "DataFrame")
            rownames(parts) <- parts[["participant_id"]]
            ## replace munged COAD with TCGA
            parts[["participant_id"]] <-
                gsub("^[A-Z]{4}", "TCGA", parts[["participant_id"]])
            coldata <- merge(
                parts, coldata, by = "participant_id",
            )
            rownames(coldata) <- coldata[["participant_id"]]
        }
        coldata <- methods::as(coldata, "DataFrame")
    }
    S4Vectors::metadata(coldata)[["columnName"]] <- columnName
    coldata
}
