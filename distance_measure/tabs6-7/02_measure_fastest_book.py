import pandas as pd

df = pd.read_csv("measure_book_week_distance", sep = "\t")
df = df.groupby(by = ["measure", "book"])["distance"].mean().reset_index()
books = set(df["book"])

src_folder = "" # This is the folder where you stored the Anobii book data (anobii-books-*.dat files)
book_metadata = pd.concat([
   pd.read_csv("%s/anobii-books-%s.dat" % (src_folder, i), sep = "\t", header = None, names = ("book", "isbn", "title", "author"), encoding = "latin1") for i in range(1, 6)
])

book_metadata = book_metadata[book_metadata["book"].isin(books)].drop_duplicates(subset = ("book",))
df = df.merge(book_metadata, on = "book")

print("Fastest")
for measure in ("lapl", "mmc", "annihil", "otp", "spls", "spla", "splc", "gft"):
   _ = df[df["measure"] == measure]
   _ = _.dropna().sort_values(by = "distance", ascending = False)
   print("%s & %s & %s\\\\" % (measure, _.iloc[0]["title"], _.iloc[0]["author"].replace("Di ", "").strip()))

print("Slowest")
for measure in ("lapl", "mmc", "annihil", "otp", "spls", "spla", "splc", "gft"):
   _ = df[df["measure"] == measure]
   _ = _.dropna().sort_values(by = "distance", ascending = True)
   print("%s & %s & %s\\\\" % (measure, _.iloc[0]["title"], _.iloc[0]["author"].replace("Di ", "").strip()))
