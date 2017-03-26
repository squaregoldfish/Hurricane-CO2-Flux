# Basic conversion functions

# Knots to metres per second
knots_to_ms <- function(knots) {
	return (as.numeric(knots) * 0.5144)
}

# Nautical miles to kilometres
nm_to_km <- function(nm) {
	return (as.numeric(nm) * 1.852)
}