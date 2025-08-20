
            library(InfIntE)
            otu_data <- read.csv("C:/MicrobiomeBot/Microbot-/temp_microbiome_9743306c.csv", row.names=1, check.names=FALSE)
            
            # Set options for Windows compatibility
            options(mc.cores = 1)
            
            # Handle Python3 dependency gracefully
            tryCatch({
                # Try to run analysis with error handling
                result <- infinte(otu_tb=otu_data, exclusion=TRUE, ncores=1, nperms=10)
            }, error = function(e) {
                # If Python3 error, try with minimal settings
                cat("Warning: Python3 dependency issue detected, using fallback mode\n")
                result <<- infinte(otu_tb=otu_data, exclusion=TRUE, ncores=1, nperms=5)
            })
            
            interactions <- result$selected_interactions
            write.csv(interactions, "C:/MicrobiomeBot/Microbot-/interactions_9743306c.csv", row.names=FALSE)
            
            # Also save enhanced interactions with metadata
            interactions_enhanced <- data.frame(
                Source = interactions$sp1,
                Target = interactions$sp2,
                Interaction_Type = interactions$lnk,
                Confidence = interactions$comp,
                Dataset_Size = "17x15",
                Optimization_Applied = "No",
                Analysis_Date = Sys.Date()
            )
            write.csv(interactions_enhanced, "C:/MicrobiomeBot/Microbot-/interactions_enhanced_9743306c.csv", row.names=FALSE)
            