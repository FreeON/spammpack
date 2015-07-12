module MatrixPair
( addSubtreeNorms
, combineZeros
, height
, ifZeroReplace
, matrixListToPair
, isZero
, MatrixPair
, MatrixTree(..)
, norm
, Norm
, pairToMatrixList
, readPairFromMatrixMarket
, setNorm
, Value
, valueNorm
, width
, writePairToMatrixMarket
) where

-- a recursive data type for matrices that efficiently encodes sparsity

import MatrixMarket (MatrixList, readFromMatrixMarket, writeToMatrixMarket)

type Value = Double ; type Norm = Double

data MatrixTree = Zero  {heightVal :: Int, widthVal :: Int} |
                  Value {normVal :: Norm, value :: Value} |
                  Rect  {heightVal :: Int, widthVal :: Int,
                         normVal :: Norm, tltree :: MatrixTree, trtree :: MatrixTree,
                                          bltree :: MatrixTree, brtree :: MatrixTree}
                  deriving (Eq, Show)

type MatrixPair = (Int, Int, MatrixTree)

height :: MatrixTree -> Int
height (Value _ _) = 1
height tree        = heightVal tree

width :: MatrixTree -> Int
width (Value _ _) = 1
width tree        = widthVal tree

norm :: MatrixTree -> Norm
norm (Zero _ _) = 0
norm tree       = normVal tree

-- setting norms

setNorm :: MatrixTree -> Norm
setNorm (Zero _ _)  = 0
setNorm (Value _ x) = valueNorm x
setNorm tree        = addSubtreeNorms . fmap setNorm $ subTrees tree

valueNorm :: Value -> Norm
valueNorm = abs

addSubtreeNorms :: [Norm] -> Norm
addSubtreeNorms = sqrt . sum . fmap (^2)

-- reading from/writing to MatrixMarket files

readPairFromMatrixMarket :: FilePath -> IO MatrixPair
readPairFromMatrixMarket filePath = readFromMatrixMarket filePath >>=
                                    (return . matrixListToPair)

writePairToMatrixMarket :: MatrixPair -> String -> FilePath -> IO ()
writePairToMatrixMarket pair format filePath =
                        writeToMatrixMarket (pairToMatrixList pair) format filePath

matrixListToPair :: MatrixList -> MatrixPair
matrixListToPair (m, n, ijxs) = (m, n, foldr addValueToTree (Zero p p) ijxs)
                                where p = nextPowOf2 $ max m n

addValueToTree :: (Int, Int, Value) -> MatrixTree -> MatrixTree
addValueToTree (i, j, x) tree = if (i, j) `inTree` tree
                                then addVal (i, j, x) tree else tree

inTree :: (Int, Int) -> MatrixTree -> Bool
inTree _ (Value _ _) = True
inTree (i, j) tree   = i <= height tree && j <= width tree

addVal :: (Int, Int, Value) -> MatrixTree -> MatrixTree

addVal (_, _, x) (Value _ _) = if x == 0 then Zero 1 1 else Value (valueNorm x) x

addVal (i, j, x) tree@(Zero h w)
       | x == 0               = tree
       | h == 1 && w == 1     = Value (valueNorm x) x
       | (i,j)   `inTree` ztl = Rect h w (valueNorm x) (addVal (i,  j,  x) ztl) ztr zbl zbr
       | (i,jr)  `inTree` ztr = Rect h w (valueNorm x) ztl (addVal (i,  jr, x) ztr) zbl zbr
       | (ib,j)  `inTree` zbl = Rect h w (valueNorm x) ztl ztr (addVal (ib,  j, x) zbl) zbr
       | (ib,jr) `inTree` zbr = Rect h w (valueNorm x) ztl ztr zbl (addVal (ib, jr, x) zbr)
       where halfh = h `div` 2 ; halfw = w `div` 2
             ib = i - halfh ; jr = j - halfw
             ztl = Zero halfh halfw
             ztr = Zero halfh (w - halfw)
             zbl = Zero (h - halfh) halfw
             zbr = Zero (h - halfh) (w - halfw)

addVal (i, j, x) (Rect h w _ tl tr bl br) = if x == 0 then ifZeroReplace newTree else newTree
       where newTree = Rect h w y newtl newtr newbl newbr
             [newtl, newtr, newbl, newbr]
                     | (i,j)   `inTree` tl = [addVal (i,  j,  x) tl, tr, bl, br]
                     | (i,jr)  `inTree` tr = [tl, addVal (i,  jr, x) tr, bl, br]
                     | (ib,j)  `inTree` bl = [tl, tr, addVal (ib,  j, x) bl, br]
                     | (ib,jr) `inTree` br = [tl, tr, bl, addVal (ib, jr, x) br]
             y = addSubtreeNorms . fmap norm $ [newtl, newtr, newbl, newbr]
             ib = i - halfh ; jr = j - halfw
             halfh = h `div` 2 ; halfw = w `div` 2

pairToMatrixList :: MatrixPair -> MatrixList
pairToMatrixList (h, w, tree) = (h, w, treeToList tree)

treeToList :: MatrixTree -> [(Int, Int, Value)]
treeToList (Zero _ _)               = []
treeToList (Value _ x)              = [(1, 1, x)]
treeToList (Rect h w _ tl tr bl br) = concat [tlijxs, fmap wshift trijxs,
                                              fmap hshift blijxs,
                                              fmap (hshift . wshift) brijxs]
     where  [tlijxs, trijxs, blijxs, brijxs] = fmap treeToList [tl, tr, bl, br]
            hshift (i, j, x) = (i + halfh, j, x)
            wshift (i, j, x) = (i, j + halfw, x)
            halfh = h `div` 2 ; halfw = w `div` 2

-- utility functions

isZero :: MatrixTree -> Bool
isZero (Zero _ _) = True
isZero _          = False

ifZeroReplace :: MatrixTree -> MatrixTree
ifZeroReplace tree@(Zero _ _)  = tree
ifZeroReplace tree@(Value _ x) = if x == 0 then Zero 1 1 else tree
ifZeroReplace tree@(Rect h w _ _ _ _ _)
                               = if all isZero (subTrees tree) then Zero h w else tree

combineZeros :: MatrixTree -> MatrixTree
combineZeros tree@(Rect h w _ _ _ _ _) = ifZeroReplace newTree
             where newTree = Rect h w x newtl newtr newbl newbr
                   [newtl, newtr, newbl, newbr] = fmap combineZeros $ subTrees tree
                   x = addSubtreeNorms . fmap norm $ [newtl, newtr, newbl, newbr]
combineZeros tree = ifZeroReplace tree

subTrees :: MatrixTree -> [MatrixTree]
subTrees (Rect _ _ _ tl tr bl br) = [tl, tr, bl, br]
subTrees _                      = []

nextPowOf2 :: Integral a => a -> a
nextPowOf2 n = head . dropWhile (< n) $ map (2^) [0..]
