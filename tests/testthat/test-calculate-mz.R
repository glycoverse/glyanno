expect_mass <- function(object, expected, ...) {
  expect_equal(round(object, 2), expected)
}

# ===== Test input types =====
test_that("calculate_mz works for glycan compositions", {
  comp <- glyrepr::glycan_composition(c(Hex = 5, HexNAc = 4, dHex = 1, NeuAc = 2))
  expect_mass(calculate_mz(comp, charge = 0), 2368.84)
})

test_that("calculate_mz works for glycan compositions as characters", {
  comp <- "Hex(5)HexNAc(4)dHex(1)NeuAc(2)"
  expect_mass(calculate_mz(comp, charge = 0), 2368.84)
})

test_that("calculate_mz works for glycan structures", {
  struc <- glyparse::parse_wurcs("WURCS=2.0/6,12,11/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O][a1221m-1a_1-5]/1-1-2-3-1-4-5-3-1-4-5-6/a4-b1_a6-l1_b4-c1_c3-d1_c6-h1_d2-e1_e4-f1_f3-g2_h2-i1_i4-j1_j3-k2")
  expect_mass(calculate_mz(struc, charge = 0), 2368.84)
})

test_that("calculate_mz works for glycan structures as characters", {
  struc <- "WURCS=2.0/6,12,11/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O][a1221m-1a_1-5]/1-1-2-3-1-4-5-3-1-4-5-6/a4-b1_a6-l1_b4-c1_c3-d1_c6-h1_d2-e1_e4-f1_f3-g2_h2-i1_i4-j1_j3-k2"
  expect_mass(calculate_mz(struc, charge = 0), 2368.84)
})

test_that("calculate_mz works for concrete monosaccharides", {
  comp <- glyrepr::glycan_composition(c(Man = 3, Gal = 2, GlcNAc = 4, Fuc = 1, Neu5Ac = 2))
  expect_mass(calculate_mz(comp, charge = 0), 2368.84)
})

# ===== Test vectorization =====
test_that("calculate_mz works for vector input", {
  comps <- c("Hex(5)HexNAc(4)dHex(1)NeuAc(2)", "HexNAc(2)Hex(3)")
  expect_mass(calculate_mz(comps, charge = 0), c(2368.84, 910.33))
})

# ===== Test derivization and mass type =====
test_that("calculate_mz works for permethylated glycans", {
  mz <- calculate_mz("HexNAc(2)Hex(3)", charge = 0, mass_dict = glyanno_mass_dict(deriv = "permethyl"))
  expect_mass(mz, 1148.59)
})

test_that("calculate_mz works for peracetylated glycans", {
  mz <- calculate_mz("HexNAc(2)Hex(3)", charge = 0, mass_dict = glyanno_mass_dict(deriv = "peracetyl"))
  expect_mass(mz, 1540.49)
})

test_that("calculate_mz works for average mass", {
  mz <- calculate_mz("HexNAc(2)Hex(3)", charge = 0, mass_dict = glyanno_mass_dict(mass_type = "average"))
  expect_mass(mz, 910.83)
})

# ===== Test charge and adduct =====
test_that("calculate_mz works for charge 1 and adduct H+", {
  mz <- calculate_mz("HexNAc(2)Hex(3)", charge = 1, adduct = "H+")
  expect_mass(mz, 911.34)
})

test_that("calculate_mz works for charge 2 and adduct H+", {
  mz <- calculate_mz("HexNAc(2)Hex(3)", charge = 2, adduct = "H+")
  expect_mass(mz, 456.17)
})

test_that("calculate_mz works for charge 1 and adduct Na+", {
  mz <- calculate_mz("HexNAc(2)Hex(3)", charge = 1, adduct = "Na+")
  expect_mass(mz, 933.32)
})

test_that("calculate_mz works for charge -1 and adduct Cl-", {
  mz <- calculate_mz("HexNAc(2)Hex(3)", charge = -1, adduct = "Cl-")
  expect_mass(mz, 945.30)
})
