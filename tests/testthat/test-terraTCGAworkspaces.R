test_that("terraTCGAworkspace works", {
    skip_if_not(
        GCPtools::gcloud_exists() &&
        identical(AnVILBase::avplatform_namespace(), "AnVILGCP") &&
        nzchar(AnVILGCP::avworkspace_name())
    )
    
    ## Blank not in validation table
    withr::with_options(
        list("terraTCGAdata.workspace" = ""),
        expect_error(
            terraTCGAworkspace()
        )
    )

    withr::with_options(
        list("terraTCGAdata.workspace" = NULL),
        expect_warning(
            terraTCGAworkspace()
        )
    )

    withr::with_options(
        list("terraTCGAdata.workspace" = "TCGA_ACC_OpenAccess_V1-0_DATA"),
        expect_identical(
            terraTCGAworkspace(),
            "TCGA_ACC_OpenAccess_V1-0_DATA"
        )
    )

    ## Replacement when input provided
    expect_warning(
        withr::with_options(
            list("terraTCGAdata.workspace" = "TCGA_ACC_OpenAccess_V1-0_DATA"),
            expect_identical(
                terraTCGAworkspace("TCGA_TGCT_OpenAccess_V1-0_DATA"),
                "TCGA_TGCT_OpenAccess_V1-0_DATA"
            )
        )
    )

    expect_error(
        terraTCGAworkspace("ABDC")
    )
})
