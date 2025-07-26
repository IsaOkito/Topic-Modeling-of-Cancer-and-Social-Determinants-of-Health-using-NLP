
install.packages(c("revtools", "tm", "topicmodels", "tidytext", "dplyr"))
library(revtools)

getwd()


# Load the .nbib file
data <- read_bibliography("pubmed-cancerANDs-set (1).nbib")

# View abstracts
head(data$abstract)

#to view all abstracts
all_abstracts <- data$abstract

cat(all_abstracts[1:3], sep = "\n\n")

# Check how many abstracts:
length(all_abstracts) 


library(revtools)       # For reading PubMed .nbib files
library(tm)             # For text preprocessing
library(topicmodels)    # For LDA topic modeling
library(tidytext)       # For tidy format text analysis
library(dplyr) 


# Convert abstracts into a text corpus (collection of text documents)
abstracts_corpus <- VCorpus(VectorSource(all_abstracts))


# Convert text to lowercase for consistency
abstracts_corpus <- tm_map(abstracts_corpus, content_transformer(tolower))

# Remove punctuation
abstracts_corpus <- tm_map(abstracts_corpus, removePunctuation)


# Remove common English stopwords (like "and", "the", "is")# TRY NOT TO REMOVE 
abstracts_corpus <- tm_map(abstracts_corpus, removeWords, stopwords("english"))

# Remove extra white spaces
abstracts_corpus <- tm_map(abstracts_corpus, stripWhitespace)

# Create a Document-Term Matrix (rows = abstracts, columns = terms/words)
dtm <- DocumentTermMatrix(abstracts_corpus)


# OPTIONAL: Remove sparse terms (words appearing in very few abstracts)
# Keeps terms that appear in at least 1% of documents (adjustable)
dtm <- removeSparseTerms(dtm, 0.999) # to SKIP THIS PART

# Set desired number of topics (change as needed)
num_topics <- 3

# Remove empty documents (rows with zero total words)
row_totals <- slam::row_sums(dtm)  # Sum of words in each document
dtm <- dtm[row_totals > 0, ]       # Keep only documents with at least 1 word


# Fit the Latent Dirichlet Allocation (LDA) model to the DTM
lda_model <- LDA(dtm, k = num_topics, control = list(seed = 1234))# NOT SET THE SEED SO TO SKIP
lda_model <- LDA(dtm, k = num_topics)

# Extract topic-word probabilities (beta matrix) from LDA model
topics <- tidy(lda_model, matrix = "beta")


# Get top 5 words per topic (words with highest probability in each topic)


top_terms <- topics %>% 
  group_by(topic) %>%  # Group by each topic

  top_n(5, beta) %>% # Select top 10 words per topic

  arrange(topic, -beta)


# Print top terms to view
print(top_terms)
print(top_terms, n = 50)   # To display all  terms (5 per topic × 10 topics)

 


write.csv(top_terms, "LDA_Topics_Terms2.csv", row.names = FALSE)

# Get the topic distribution for each document
doc_topics <- tidy(lda_model, matrix = "gamma")

# View or export document-topic probabilities
head(doc_topics)

# Optional: Save document-topic assignment
write.csv(doc_topics, "LDA_Document_Topic_Distribution2.csv", row.names = FALSE)


#visualization

library(ggplot2)

my_plot <- ggplot(top_terms, aes(x = reorder(term, beta), y = beta, fill = factor(topic))) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free", ncol = 2) +  # Change ncol as desired
  coord_flip() +
  labs(title = "Top 5 Words per Topic",
       x = "Top Words",
       y = "Word Importance (Beta Value)") +
  theme_minimal()

# Display plot
print(my_plot)



# Save as PNG
ggsave("topic_model_plot.png", plot = my_plot, width = 10, height = 8, dpi = 300)




