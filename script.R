
            library(InfIntE)
            otu_data <- read.csv("C:/MicrobiomeBot/Microbot-/input.csv", row.names=1, check.names=FALSE)
            
            # Set options for Windows compatibility and speed optimization
            options(mc.cores = 1)
            
            # Handle Python3 dependency gracefully with optimized parameters
            tryCatch({
                # Try to run analysis with error handling
                result <- infinte(otu_tb=otu_data, exclusion=TRUE, ncores=1, nperms=5, lambda.min.ratio=0.1)
            }, error = function(e) {
                # If Python3 error, try with minimal settings
                cat("Warning: Python3 dependency issue detected, using fallback mode\n")
                result <<- infinte(otu_tb=otu_data, exclusion=TRUE, ncores=1, nperms=3, lambda.min.ratio=0.1)
            })
            
            interactions <- result$selected_interactions
            write.csv(interactions, "C:/MicrobiomeBot/Microbot-/interactions.csv", row.names=FALSE)
            
            # Also save enhanced interactions with metadata
            interactions_enhanced <- data.frame(
                Species1 = interactions$sp1,
                Species2 = interactions$sp2,
                Interaction_Type = interactions$lnk,
                Compression_Value = interactions$comp,
                Interaction_Strength = cut(interactions$comp, 
                                          breaks = quantile(interactions$comp, c(0, 0.33, 0.66, 1)), 
                                          labels = c("Weak", "Medium", "Strong"),
                                          include.lowest = TRUE),
                Direction = ifelse(grepl("up|app", interactions$lnk), "Positive", "Negative"),
                Timestamp = Sys.time()
            )
            
            write.csv(interactions_enhanced, "C:/MicrobiomeBot/Microbot-/interactions_enhanced.csv", row.names=FALSE)
            