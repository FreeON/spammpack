module MatrixTree
( addSubtreeNorms
, combineZeros
, height
, ifZeroReplace
, indexedListToTree
, isZero
, MatrixTree(..)
, norm
, Norm
, readTreeFromMatrixMarket
, rectOrder
, setNorm
, treeToIndexedList
, Value
, valueNorm
, width
, writeTreeToMatrixMarket
) where

-- a recursive data type for matrices that efficiently encodes sparsity

import Data.List (partition)
import MatrixMarket (IndexedList, readFromMatrixMarket, writeToMatrixMarket)

type Value = Double ; type Norm = Double

data MatrixTree = Zero  {top :: Int, left :: Int, heightVal :: Int, widthVal :: Int} |
                  Value {top :: Int, left :: Int,
                         normVal :: Norm, value :: Value} |
                  Rect  {top :: Int, left :: Int, heightVal :: Int, widthVal :: Int,
                         normVal :: Norm, tltree :: MatrixTree, trtree :: MatrixTree,
                                          bltree :: MatrixTree, brtree :: MatrixTree}
                  deriving (Eq, Show)

height :: MatrixTree -> Int
height (Value _ _ _ _) = 1
height tree            = heightVal tree

width :: MatrixTree -> Int
width (Value _ _ _ _) = 1
width tree            = widthVal tree

norm :: MatrixTree -> Norm
norm (Zero _ _ _ _) = 0
norm tree           = normVal tree

-- a Zero with height and width > 0 is a block with all zero entries;
-- a Zero with height or width = 0 is a block with no entries at all

-- setting norms

setNorm :: MatrixTree -> Norm
setNorm (Zero _ _ _ _)               = 0
setNorm (Value _ _ _ x)              = valueNorm x
setNorm (Rect _ _ _ _ _ tl tr bl br) = addSubtreeNorms . fmap setNorm $ [tl, tr, bl, br]

valueNorm :: Value -> Norm
valueNorm = abs

addSubtreeNorms :: [Norm] -> Norm
addSubtreeNorms = sqrt . sum . fmap (^2)

-- reading from/writing to MatrixMarket files

readTreeFromMatrixMarket :: FilePath -> IO MatrixTree
readTreeFromMatrixMarket filePath = readFromMatrixMarket filePath >>=
                                    (return . indexedListToTree)

writeTreeToMatrixMarket :: MatrixTree -> String -> FilePath -> IO ()
writeTreeToMatrixMarket tree format filePath =
                        writeToMatrixMarket (treeToIndexedList tree) format filePath

indexedListToTree :: IndexedList -> MatrixTree
indexedListToTree (m, n, ijxs) = foldr addValueToTree (Zero 1 1 m n) ijxs

addValueToTree :: (Int, Int, Value) -> MatrixTree -> MatrixTree
addValueToTree (i, j, x) tree = if (i, j) `inTree` tree then addVal (i, j, x) tree else tree

inTree :: (Int, Int) -> MatrixTree -> Bool
inTree (i, j) (Value m n _ _) = i == m && j == n
inTree (i, j) tree            = inRange (i, j) (top tree, left tree, height tree, width tree)

inRange :: (Int, Int) -> (Int, Int, Int, Int) -> Bool
inRange (i, j) (t, l, h, w) = i >= t && i <= t + h - 1 && j >= l && j <= l + w - 1

addVal :: (Int, Int, Value) -> MatrixTree -> MatrixTree

addVal (i, j, x) (Value _ _ _ _) = if x == 0 then Zero i j 1 1 else Value i j (valueNorm x) x

addVal (i, j, x) tree@(Zero t l h w)
       | x == 0             = tree
       | [h,w] == [1,1]     = Value t l (valueNorm x) x
       | (i,j) `inTree` ztl = Rect t l h w (valueNorm x) (addVal (i, j, x) ztl) ztr zbl zbr
       | (i,j) `inTree` ztr = Rect t l h w (valueNorm x) ztl (addVal (i, j, x) ztr) zbl zbr
       | (i,j) `inTree` zbl = Rect t l h w (valueNorm x) ztl ztr (addVal (i, j, x) zbl) zbr
       | (i,j) `inTree` zbr = Rect t l h w (valueNorm x) ztl ztr zbl (addVal (i, j, x) zbr)
       where ztl = Zero t l halfh halfw
             ztr = Zero t (l + halfw) halfh (w - halfw)
             zbl = Zero (t + halfh) l (h - halfh) halfw
             zbr = Zero (t + halfh) (l + halfw) (h - halfh) (w - halfw)
             halfh = h `div` 2 ; halfw = w `div` 2

addVal (i, j, x) (Rect t l h w _ tl tr bl br) = if x == 0 then ifZeroReplace newTree else newTree
       where newTree = Rect t l h w y newtl newtr newbl newbr
             [newtl, newtr, newbl, newbr]
                     | (i,j) `inTree` tl = [addVal (i, j, x) tl, tr, bl, br]
                     | (i,j) `inTree` tr = [tl, addVal (i, j, x) tr, bl, br]
                     | (i,j) `inTree` bl = [tl, tr, addVal (i, j, x) bl, br]
                     | (i,j) `inTree` br = [tl, tr, bl, addVal (i, j, x) br]
             y = addSubtreeNorms . fmap norm $ [newtl, newtr, newbl, newbr]

treeToIndexedList :: MatrixTree -> IndexedList
treeToIndexedList (Zero _ _ h w)               = (h, w, [])
treeToIndexedList (Value _ _ _ x)              = (1, 1, [(1, 1, x)])
treeToIndexedList (Rect _ _ h w _ tl tr bl br) = (h, w, ijxs)
      where ijxs = concat [tlijxs, fmap wshift trijxs,
                           fmap hshift blijxs, fmap (hshift . wshift) brijxs]
            [tlijxs, trijxs, blijxs, brijxs] = fmap (third . treeToIndexedList) [tl, tr, bl, br]
            third (a, b, c) = c
            hshift (i, j, x) = (i + halfh, j, x)
            wshift (i, j, x) = (i, j + halfw, x)
            halfh = h `div` 2 ; halfw = w `div` 2

-- utility functions

isZero :: MatrixTree -> Bool
isZero (Zero _ _ _ _) = True
isZero _              = False

ifZeroReplace :: MatrixTree -> MatrixTree
ifZeroReplace tree@(Zero _ _ _ _) = tree
ifZeroReplace tree@(Value i j _ x) = if x == 0 then Zero i j 1 1 else tree
ifZeroReplace tree@(Rect t l h w _ tl tr bl br)
              | null nonZeroTrees       = Zero t l h w
              | length nonZeroTrees > 1 = tree
              | otherwise               = if alone then head nonZeroTrees else tree
              where (zeroTrees, nonZeroTrees) = partition isZero [tl, tr, bl, br]
                    alone = any (\tr -> height tr == 0 && width tr == 0) zeroTrees

combineZeros :: MatrixTree -> MatrixTree
combineZeros (Rect t l h w _ tl tr bl br) =
             ifZeroReplace (Rect t l h w x newtl newtr newbl newbr)
             where [newtl, newtr, newbl, newbr] = fmap combineZeros [tl, tr, bl, br]
                   x = addSubtreeNorms . fmap norm $ [newtl, newtr, newbl, newbr]
combineZeros tree = ifZeroReplace tree

rectOrder :: MatrixTree -> MatrixTree
rectOrder tree@(Zero _ _ _ _) = tree
rectOrder tree@(Value _ _ _ _) = tree
rectOrder tree@(Rect t l h w x tl tr bl br)
          | fmap height [bl, br] == [0, 0] = Rect t l h w x bl br tl tr
          | fmap width  [tr, br] == [0, 0] = Rect t l h w x tr tl br bl
          | otherwise                      = tree
