library(chatgpt)
# Key

key1 <- "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"

Sys.setenv(OPENAI_API_KEY = key1)

cat(ask_chatgpt("What do you think about the R language?"))
