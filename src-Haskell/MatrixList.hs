module MatrixList
( MatrixList
, mEntryCol
, mEntryRow
, mEntryVal
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

mEntryRow :: MListEntry a -> Int
mEntryRow (i, _, _) = i

mEntryCol :: MListEntry a -> Int
mEntryCol (_, j, _) = j

mEntryVal :: MListEntry a -> a
mEntryVal (_, _, x) = x
