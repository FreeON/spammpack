module VectorList
( entryIndex
, entryVal
, VectorList
, VListEntry
, vListEntries
, vListLength
) where

type VListEntry a = (Int, a)
type VectorList a = (Int, [VListEntry a])

-- vector length and list of entries (index, value) ;
-- including zero-value entries is optional

vListLength :: VectorList a -> Int
vListLength (l, _) = l

vListEntries :: VectorList a -> [VListEntry a]
vListEntries (_, ixs) = ixs

entryIndex :: VListEntry a -> Int
entryIndex (i, _) = i

entryVal :: VListEntry a -> a
entryVal (_, x) = x
