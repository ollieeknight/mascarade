library(ggplot2)

make_mask_table <- function() {
    data("exampleMascarade")
    generateMask(dims = exampleMascarade$dims,
                 clusters = exampleMascarade$clusters,
                 gridSize = 50)
}

get_layer_colors <- function(res) {
    # The geom layer is the second element (first is coord_cartesian)
    # Note: "colour" is the ggplot2 internal field name
    res[[2]]$aes_params$colour
}

test_that("fancyMask uses default hue palette when cols is NULL", {
    mt <- make_mask_table()
    res <- fancyMask(mt)
    colors <- get_layer_colors(res)
    clusterLevels <- levels(mt$cluster)
    expected_pal <- setNames(scales::hue_pal()(length(clusterLevels)), clusterLevels)
    expect_equal(unname(colors), unname(expected_pal[mt$cluster]))
})

test_that("fancyMask applies single color to all clusters", {
    mt <- make_mask_table()
    res <- fancyMask(mt, cols = "red")
    colors <- get_layer_colors(res)
    expect_true(all(colors == "red"))
})

test_that("fancyMask maps unnamed color vector by factor level order", {
    mt <- make_mask_table()
    clusterLevels <- levels(mt$cluster)
    col_vec <- rep("grey50", length(clusterLevels))
    col_vec[1] <- "blue"
    res <- fancyMask(mt, cols = col_vec)
    colors <- get_layer_colors(res)
    expected_pal <- setNames(col_vec, clusterLevels)
    expect_equal(unname(colors), unname(expected_pal[mt$cluster]))
})

test_that("fancyMask maps named color vector by cluster name", {
    mt <- make_mask_table()
    clusterLevels <- levels(mt$cluster)
    # Provide in reverse order to verify matching is by name, not position
    col_vec <- setNames(scales::hue_pal()(length(clusterLevels)),
                        rev(clusterLevels))
    res <- fancyMask(mt, cols = col_vec)
    colors <- get_layer_colors(res)
    expect_equal(unname(colors), unname(col_vec[as.character(mt$cluster)]))
})

test_that("fancyMask errors on wrong-length unnamed cols", {
    mt <- make_mask_table()
    expect_error(fancyMask(mt, cols = c("red", "blue")),
                 "must equal the number of clusters")
})

test_that("fancyMask errors on named cols missing a cluster", {
    mt <- make_mask_table()
    clusterLevels <- levels(mt$cluster)
    # Provide only the first two clusters, omitting the rest
    col_vec <- setNames(c("red", "blue"), clusterLevels[1:2])
    expect_error(fancyMask(mt, cols = col_vec),
                 "missing entries for cluster")
})

test_that("fancyMask renders without error with custom cols", {
    mt <- make_mask_table()
    data("exampleMascarade")
    p <- ggplot(do.call(cbind, exampleMascarade)) +
        geom_point(aes(x = UMAP_1, y = UMAP_2)) +
        fancyMask(mt, cols = "red") +
        theme_classic()
    pf <- tempfile(fileext = ".pdf")
    expect_no_error(ggsave(p, file = pf, width = 5, height = 4))
    expect_true(file.exists(pf))
})
