# Hypothetical data
df <- data.frame(
  Sample_ID = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3"),
  Collection_Date = c(
    "01-06-2022", "02-09-2022", "01-01-2023",
    "15-05-2022", "15-08-2022",
    "10-07-2022", "10-10-2022", "10-01-2023", "10-04-2023"
  ),
  stringsAsFactors = FALSE
)

# Convert to Date type
df$Collection_Date <- as.Date(df$Collection_Date, format = "%d-%m-%Y")

# Functions to convert date to integer depending on frequency -----------------

# Annual
annual_time <- function(date) {
  as.integer(format(date, "%Y"))
}

# Semi-annual
semiannual_time <- function(date) {
  yr <- as.integer(format(date, "%Y"))
  half <- ifelse(as.integer(format(date, "%m")) <= 6, 0, 0.5)
  as.integer(( (yr + half) * 2 ) + 1)
}

# Quarterly
quarterly_time <- function(date) {
  yr <- as.integer(format(date, "%Y"))
  qtr <- ((as.integer(format(date, "%m")) - 1) %/% 3) * 0.25
  as.integer(((yr + qtr) * 4) + 1)
}

# Monthly
monthly_time <- function(date) {
  yr <- as.integer(format(date, "%Y"))
  mo <- (as.integer(format(date, "%m")) - 1) / 12
  as.integer(((yr + mo) * 12) + 1)
}

# ---------------------------------------------------------------------------

# Apply conversions
df$Annual_Int     <- annual_time(df$Collection_Date)
df$SemiAnnual_Int <- semiannual_time(df$Collection_Date)
df$Quarterly_Int  <- quarterly_time(df$Collection_Date)
df$Monthly_Int    <- monthly_time(df$Collection_Date)

print(df)
