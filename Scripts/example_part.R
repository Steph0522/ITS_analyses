doubs.spe <- read.csv ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/DoubsSpe.csv', row.names = 1)
doubs.env <- read.csv ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/DoubsEnv.csv', row.names = 1)
doubs.spa <- read.csv ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/DoubsSpa.csv', row.names = 1)

spe <- doubs.spe

env.z <- decostand(doubs.env, method = "standardize")
env.z <- subset(env.z, select = -das)
spe.hel <- decostand(spe, method = "hellinger")
spe.rda <- rda(spe.hel ~ ., data = env.z)

# Model the effect of all environmental variables on fish
# community composition
spe.rda <- rda(spe.hel ~ ., data = env.z)
summary(spe.rda)

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(spe.hel ~ 1, data = env.z), # lower model limit (simple!)
                      scope = formula(spe.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) # change to TRUE to see the selection process!
# Check the new model with forward-selected variables
fwd.sel$call

# Write our new model
spe.rda.signif <- rda(spe.hel ~ alt + oxy + dbo, data = env.z)
# check the adjusted R2 (corrected for the number of
# explanatory variables)
RsquareAdj(spe.rda.signif)

# Subset environmental data into topography variables and
# chemistry variables
env.topo <- subset(env.z, select = c(alt, pen, deb))
env.chem <- subset(env.z, select = c(pH, dur, pho, nit, amm,
                                     oxy, dbo))

# Run a partial RDA
spe.partial.rda <- rda(spe.hel, env.chem, env.topo)

# Partition the variation in fish community composition
spe.part.all <- varpart(spe.hel, env.chem, env.topo)
spe.part.all$part  # access results!

# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("Chem", "Topo"), # name the partitions
     bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
