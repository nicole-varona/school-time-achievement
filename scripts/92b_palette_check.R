# ==============================================================================
# 92b_palette_check.R — Colour-vision validation for the figure palette
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Outputs : a table on stdout. Writes nothing, and is NOT part of run_all.R: it is
#           a design-time check, run by hand whenever the palette in 90_style.R
#           changes, so that a colour choice is never made by eye.
#
# The floor the project uses is a separation of 8 (OKLab distance x100) under
# normal vision, deuteranopia and protanopia.
# ==============================================================================

# Separación percibida entre dos colores bajo visión normal, deuteranopía y
# protanopía. Simulación de Machado, Oliveira y Fernandes (2009) sobre RGB lineal,
# distancia euclídea en OKLab x100. El piso que fijó el proyecto es 8.
srgb_to_lin <- function(c) ifelse(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055)^2.4)

M_DEUT <- matrix(c( 0.367322, 0.860646, -0.227968,
                    0.280085, 0.672501,  0.047413,
                   -0.011820, 0.042940,  0.968881), 3, 3, byrow = TRUE)
M_PROT <- matrix(c( 0.152286, 1.052583, -0.204868,
                    0.114503, 0.786281,  0.099216,
                   -0.003882,-0.048116,  1.051998), 3, 3, byrow = TRUE)

simulate <- function(hex, M = NULL) {
  v <- srgb_to_lin(as.numeric(grDevices::col2rgb(hex)) / 255)
  if (!is.null(M)) v <- as.numeric(M %*% v)
  pmin(pmax(v, 0), 1)
}

lin_to_oklab <- function(v) {
  l <- 0.4122214708*v[1] + 0.5363325363*v[2] + 0.0514459929*v[3]
  m <- 0.2119034982*v[1] + 0.6806995451*v[2] + 0.1073969566*v[3]
  s <- 0.0883024619*v[1] + 0.2817188376*v[2] + 0.6299787005*v[3]
  l_ <- l^(1/3); m_ <- m^(1/3); s_ <- s^(1/3)
  c(0.2104542553*l_ + 0.7936177850*m_ - 0.0040720468*s_,
    1.9779984951*l_ - 2.4285922050*m_ + 0.4505937099*s_,
    0.0259040371*l_ + 0.7827717662*m_ - 0.8086757660*s_)
}

sep <- function(a, b, M = NULL)
  round(100 * sqrt(sum((lin_to_oklab(simulate(a, M)) - lin_to_oklab(simulate(b, M)))^2)), 1)

check <- function(name, a, b) {
  s <- c(normal = sep(a, b), deuteranopia = sep(a, b, M_DEUT), protanopia = sep(a, b, M_PROT))
  cat(sprintf("%-34s %-9s %-9s  normal %5.1f  deut %5.1f  prot %5.1f   %s\n",
              name, a, b, s[1], s[2], s[3],
              if (min(s) >= 8) "PASA" else "FALLA"))
  invisible(s)
}

cat("--- referencia: el par actual ---\n")
check("Okabe-Ito azul / naranja", "#0072B2", "#E69F00")

cat("\n--- candidatos en rojos/rosados, grises y celestes ---\n")
check("RdBu fuerte (rojo / azul)",      "#B2182B", "#2166AC")
check("RdBu medio (rojo / celeste)",    "#D6604D", "#4393C3")
check("RdBu claro (rosado / celeste)",  "#F4A582", "#92C5DE")
check("Rojo oscuro / celeste claro",    "#A50F15", "#9ECAE1")
check("Carmin / gris azulado",          "#C0392B", "#7F8C9A")
check("Rosado / gris oscuro",           "#E7969C", "#636363")
check("Borgona / celeste",              "#8C2D04", "#74A9CF")
check("Rojo teja / gris pizarra",       "#CB4B41", "#6E7B8B")

cat("\n--- celeste algo mas profundo, para que la linea se lea sobre blanco ---\n")
check("Rojo oscuro / celeste medio", "#A50F15", "#6BAED6")
check("Borgona / celeste medio",     "#8C2D04", "#6BAED6")

cat("\n--- secuencia categorica descartada (se conserva por el registro) ---\n")
CAT <- c("#A50F15", "#6BAED6", "#6E7B8B", "#E7969C", "#08519C", "#BDBDBD")
nm  <- c("rojo oscuro", "celeste", "gris pizarra", "rosado", "azul profundo", "gris claro")
for (i in 1:(length(CAT) - 1))
  check(sprintf("%s / %s", nm[i], nm[i + 1]), CAT[i], CAT[i + 1])

# ==============================================================================
# LA PALETA EN USO, leida de 90_style.R y no copiada aca.
#
# Hasta el 13/08 este bloque validaba una secuencia tipeada a mano que ya no era
# la del proyecto: difería de PAL_CAT en dos de seis colores (#A50F15 por
# #8C2D04 y #6E7B8B por #636363). O sea que correr el chequeo "validaba" una
# paleta que ningun grafico usaba, que es peor que no correrlo. Ahora se sourcea
# 90_style.R, de modo que el chequeo no puede volver a separarse de lo validado.
# ==============================================================================
suppressPackageStartupMessages(source(here::here("scripts", "90_style.R")))

cat("\n=== PALETA EN USO: el par de materias (PAL) ===\n")
check("matematica / lengua", PAL$mathematics, PAL$language)

cat("\n=== PALETA EN USO: PAL_CAT, pares adyacentes ===\n")
for (i in seq_len(length(PAL_CAT) - 1))
  check(sprintf("PAL_CAT[%d] / PAL_CAT[%d]", i, i + 1), PAL_CAT[i], PAL_CAT[i + 1])

# Los tres primeros son los que se usan en la practica (transiciones usa 1-2,
# contaminacion usa 1-2-3), asi que van todos contra todos y no solo adyacentes.
cat("\n=== PALETA EN USO: los tres primeros de PAL_CAT, todos contra todos ===\n")
for (i in 1:2) for (j in (i + 1):3)
  check(sprintf("PAL_CAT[%d] / PAL_CAT[%d]", i, j), PAL_CAT[i], PAL_CAT[j])

# El piso lo verifica el script, no el ojo de quien lee la tabla.
pairs_used <- c(list(c(PAL$mathematics, PAL$language)),
                lapply(seq_len(length(PAL_CAT) - 1), function(i) PAL_CAT[i:(i + 1)]),
                list(PAL_CAT[c(1, 3)]), list(PAL_CAT[c(2, 3)]))
worst <- min(vapply(pairs_used, function(p)
  min(sep(p[1], p[2]), sep(p[1], p[2], M_DEUT), sep(p[1], p[2], M_PROT)), numeric(1)))
cat(sprintf("\n=== VEREDICTO: separacion minima sobre los pares en uso = %.1f (piso 8) -> %s ===\n",
            worst, if (worst >= 8) "PASA" else "FALLA"))
if (worst < 8)
  stop("la paleta de 90_style.R cae por debajo del piso de separacion del proyecto")
