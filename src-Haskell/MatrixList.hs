module MatrixList
( entryCol
, entryRow
, entryVal
, MatrixList
, MListEntry
, mListEntries
, mListHeight
, mListWidth
) where

type MListEntry a = (Int, Int, a)
type MatrixList a = (Int, Int, [MListEntry a])

-- matrix height, width, and list of entries (row, column, value) ;
-- including zero-value entries is optional

mListHeight :: MatrixList a -> Int
mListHeight (h, _, _) = h

mListWidth :: MatrixList a -> Int
mListWidth (_, w, _) = w

mListEntries :: MatrixList a -> [MListEntry a]
mListEntries (_, _, ijxs) = ijxs

entryRow :: MListEntry a -> Int
entryRow (i, _, _) = i

entryCol :: MListEntry a -> Int
entryCol (_, j, _) = j

entryVal :: MListEntry a -> a
entryVal (_, _, x) = x
