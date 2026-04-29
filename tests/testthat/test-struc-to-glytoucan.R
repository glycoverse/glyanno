glytoucan_response <- function(status_code = 200, body = "{}") {
  httr2::response(
    status_code = status_code,
    headers = list("content-type" = "application/json"),
    body = charToRaw(body)
  )
}

test_that("struc_to_glytoucan maps successful responses to accessions", {
  responses <- list(
    glytoucan_response(body = '{"id":"G00001MO"}'),
    glytoucan_response(body = '{"id":"G00002MO"}')
  )
  performed_urls <- character()
  response_index <- 0
  local_mocked_bindings(
    req_perform = function(req) {
      response_index <<- response_index + 1
      performed_urls <<- c(performed_urls, req$url)
      responses[[response_index]]
    },
    .package = "httr2"
  )

  strucs <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "GalNAc(a1-"
  ))
  accessions <- suppressMessages(struc_to_glytoucan(strucs))

  expect_equal(accessions, c("G00001MO", "G00002MO"))
  expect_equal(length(performed_urls), 2)
  expect_true(all(grepl("iupaccondensed2wurcs", performed_urls, fixed = TRUE)))
})

test_that("struc_to_glytoucan returns NA for unsuccessful responses", {
  responses <- list(
    glytoucan_response(status_code = 404),
    glytoucan_response(body = '{"wurcs":"WURCS=2.0/1,1,0/[...]/1/"}')
  )
  response_index <- 0
  local_mocked_bindings(
    req_perform = function(req) {
      response_index <<- response_index + 1
      responses[[response_index]]
    },
    .package = "httr2"
  )

  strucs <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "GalNAc(a1-"
  ))
  accessions <- suppressMessages(struc_to_glytoucan(strucs))

  expect_equal(accessions, c(NA_character_, NA_character_))
})

test_that("struc_to_glytoucan accepts character structure input", {
  local_mocked_bindings(
    req_perform = function(req) {
      glytoucan_response(body = '{"id":"G00003MO"}')
    },
    .package = "httr2"
  )

  accessions <- suppressMessages(struc_to_glytoucan("GalNAc(a1-"))

  expect_equal(accessions, "G00003MO")
})

test_that("struc_to_glytoucan accepts empty structure input", {
  local_mocked_bindings(
    req_perform = function(req) {
      stop("No request should be performed for empty input.", call. = FALSE)
    },
    .package = "httr2"
  )

  accessions <- suppressMessages(struc_to_glytoucan(glyrepr::glycan_structure()))

  expect_equal(accessions, character())
})
