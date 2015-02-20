module MatrixTree
( Value
, Norm
, MatrixTree(..)
, RowList
, getNorm
, valueNorm
, addSubtreeNorms
, isZero
, ifZeroReplace
, combineZeros
, getHeight
, getWidth
, rectOrder
, treeToRowList
, rowListToTree
) where

-- a recursive data type for matrices that efficiently encodes sparsity

import Data.List (transpose)

type Value = Double ; type Norm = Double

data MatrixTree = Zero  {top :: Int, left :: Int, height :: Int, width :: Int} |
                  Value {top :: Int, left :: Int,
                         norm :: Norm, value :: Value} |
                  Rect  {top :: Int, left :: Int, height :: Int, width :: Int,
                         norm :: Norm, tltree :: MatrixTree, trtree :: MatrixTree,
                                       bltree :: MatrixTree, brtree :: MatrixTree}
                  deriving (Eq, Show)

-- a Zero with height and width > 0 is a block with all zero entries;
-- a Zero with height or width = 0 is a block with no entries at all

type RowList = [[Value]]

-- setting and accessing norms

getNorm :: MatrixTree -> Norm
getNorm (Zero _ _ _ _) = 0
getNorm tree = norm tree

valueNorm :: Value -> Norm
valueNorm = abs

addSubtreeNorms :: [Norm] -> Norm
addSubtreeNorms = sqrt . sum . fmap (^2)

-- MatrixTree utility functions

isZero :: MatrixTree -> Bool
isZero (Zero _ _ _ _) = True
isZero _              = False

ifZeroReplace :: MatrixTree -> MatrixTree
ifZeroReplace tree@(Zero _ _ _ _) = tree
ifZeroReplace tree@(Value i j _ x) = if x == 0 then Zero i j 1 1 else tree
ifZeroReplace tree@(Rect t l h w _ tl tr bl br)
              | null nonZeroTrees        = Zero t l h w
              | length nonZeroTrees == 1 = head nonZeroTrees
              | otherwise                = tree
              where nonZeroTrees = filter (not . isZero) [tl, tr, bl, br]

combineZeros :: MatrixTree -> MatrixTree
combineZeros (Rect t l h w _ tl tr bl br) =
             ifZeroReplace (Rect t l h w x newtl newtr newbl newbr)
             where [newtl, newtr, newbl, newbr] = fmap combineZeros [tl, tr, bl, br]
                   x = addSubtreeNorms . fmap getNorm $ [newtl, newtr, newbl, newbr]
combineZeros tree = ifZeroReplace tree

getHeight :: MatrixTree -> Int
getHeight (Value _ _ _ _) = 1
getHeight tree = height tree

getWidth :: MatrixTree -> Int
getWidth (Value _ _ _ _) = 1
getWidth tree = width tree

rectOrder :: MatrixTree -> MatrixTree
rectOrder tree@(Zero _ _ _ _) = tree
rectOrder tree@(Value _ _ _ _) = tree
rectOrder tree@(Rect t l h w x tl tr bl br)
          | fmap getHeight [bl, br] == [0, 0] = Rect t l h w x bl br tl tr
          | fmap getWidth [tr, br]  == [0, 0] = Rect t l h w x tr tl br bl
          | otherwise                         = tree

-- converting between data types

rowListToTree :: RowList -> MatrixTree
rowListToTree = snd . fillFromList . addListSize

treeToRowList :: MatrixTree -> RowList
treeToRowList (Zero t l h w)               = replicate h $ replicate w 0
treeToRowList (Value _ _ _ x)              = [[x]]
treeToRowList (Rect _ _ _ _ _ tl tr bl br) =
              zipWith (++) (treeToRowList tl ++ treeToRowList bl)
                           (treeToRowList tr ++ treeToRowList br)

addListSize :: RowList -> (RowList, Int, Int, Int, Int)
addListSize xss = (xss, 1, 1, length xss, length . head $ xss)

fillFromList :: (RowList, Int, Int, Int, Int) -> (Norm, MatrixTree)
fillFromList ([[x]], i, j, _, _) = if x == 0 then (0, Zero i j 1 1)
                                   else (valueNorm x, Value i j (valueNorm x) x)
fillFromList (xss, t, l, h, w)
             | all0 xss  = (0, Zero t l h w)
             | otherwise = (rectNorm, Rect t l h w rectNorm (tree tl) (tree tr)
                                                            (tree bl) (tree br))
             where all0 = and . fmap (== 0) . concat
                   [xsstl, xsstr, xssbl, xssbr] = partitionList xss
                   tl = (xsstl, t, l, h `div` 2, w `div` 2)
                   tr = (xsstr, t, l + w `div` 2, h `div` 2, w - w `div` 2)
                   bl = (xssbl, t + h `div` 2, l, h - h `div` 2, w `div` 2)
                   br = (xssbr, t + h `div` 2, l + w `div` 2, h - h `div` 2, w - w `div` 2)
                   rectNorm = calcNorm [tl, tr, bl, br]
                   calcNorm = addSubtreeNorms . fmap (fst . fillFromList)
                   tree = snd . fillFromList

partitionList :: [[a]] -> [[[a]]]
partitionList xss = fmap transpose [tlxss, trxss, blxss, brxss]
                    where (txss, bxss)   = halveList xss
                          (tlxss, trxss) = halveList . transpose $ txss
                          (blxss, brxss) = halveList . transpose $ bxss

halveList :: [a] -> ([a],[a])
halveList xs = splitAt (length xs `div` 2) xs