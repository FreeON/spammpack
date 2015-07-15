module MatrixTree
( addSubtreeNorms
, ifZeroReplace
, isZero
, matrixListToTree
, MatrixTree
, mmReadTree
, mmWriteTree
, MTree(..)
, nextPowOf2
, norm
, Norm
, setNorm
, size
, Size
, treeToMatrixList
, Value
, valueNorm
) where

-- a recursive matrix data type that efficiently encodes sparsity

import MatrixMarket (MatrixList, mmReadFile, mmWriteFile)

type Size = Int ; type Value = Double ; type Norm = Double

data MTree = Zero Size |
             Leaf Norm Value |
             Square Size Norm MTree MTree MTree MTree
             deriving (Eq, Show)

type MatrixTree = (Int, Int, MTree)

size :: MTree -> Size
size (Zero s)             = s
size (Leaf _ _)           = 1
size (Square s _ _ _ _ _) = s

norm :: MTree -> Norm
norm (Zero _)             = 0
norm (Leaf n _)           = n
norm (Square _ n _ _ _ _) = n

subTrees :: MTree -> [MTree]
subTrees (Square _ _ tl tr bl br) = [tl, tr, bl, br]
subTrees _ = []

-- setting norms

setNorm :: MTree -> Norm
setNorm (Zero _)   = 0
setNorm (Leaf _ x) = valueNorm x
setNorm tree       = addSubtreeNorms . fmap setNorm $ subTrees tree

valueNorm :: Value -> Norm
valueNorm = abs

addSubtreeNorms :: [Norm] -> Norm
addSubtreeNorms = sqrt . sum . fmap (^2)

-- reading from/writing to MatrixMarket files

mmReadTree :: FilePath -> IO MatrixTree
mmReadTree filePath = mmReadFile filePath >>= (return . matrixListToTree)

mmWriteTree :: MatrixTree -> String -> FilePath -> IO ()
mmWriteTree tree format filePath =
            mmWriteFile (treeToMatrixList tree) format filePath

matrixListToTree :: MatrixList -> MatrixTree
matrixListToTree (m, n, ijxs) = (m, n, foldr addVal (Zero p) ijxs)
                                where p = nextPowOf2 $ max m n

addVal :: (Int, Int, Value) -> MTree -> MTree
addVal (_, _, x) (Leaf _ _) = if x == 0 then Zero 1 else Leaf (valueNorm x) x
addVal (i, j, x) tree@(Zero s)
 | x == 0         = tree
 | s == 1         = Leaf (valueNorm x) x
 | within [i,j]   = Square s (valueNorm x) (addVal (i,  j,  x) zro) zro zro zro
 | within [i,jr]  = Square s (valueNorm x) zro (addVal (i,  jr, x) zro) zro zro
 | within [ib,j]  = Square s (valueNorm x) zro zro (addVal (ib,  j, x) zro) zro
 | within [ib,jr] = Square s (valueNorm x) zro zro zro (addVal (ib, jr, x) zro)
 where halfs = s `div` 2
       within = all (<= halfs)
       ib = i - halfs ; jr = j - halfs
       zro = Zero halfs
addVal (i, j, x) (Square s _ tl tr bl br) = if x == 0 then
                                            ifZeroReplace newTree else newTree
       where newTree = Square s y newtl newtr newbl newbr
             [newtl, newtr, newbl, newbr]
                     | within [i,j]   = [addVal (i,  j,  x) tl, tr, bl, br]
                     | within [i,jr]  = [tl, addVal (i,  jr, x) tr, bl, br]
                     | within [ib,j]  = [tl, tr, addVal (ib,  j, x) bl, br]
                     | within [ib,jr] = [tl, tr, bl, addVal (ib, jr, x) br]
             y = addSubtreeNorms . fmap norm $ [newtl, newtr, newbl, newbr]
             within = all (<= halfs)
             ib = i - halfs ; jr = j - halfs
             halfs = s `div` 2

treeToMatrixList :: MatrixTree -> MatrixList
treeToMatrixList (h, w, mTree) = (h, w, mTreeToList mTree)

mTreeToList :: MTree -> [(Int, Int, Value)]
mTreeToList (Zero _)                 = []
mTreeToList (Leaf _ x)               = [(1, 1, x)]
mTreeToList (Square s _ tl tr bl br) = concat [tlijxs, fmap wshift trijxs,
                                               fmap hshift blijxs,
                                               fmap (hshift . wshift) brijxs]
      where [tlijxs, trijxs, blijxs, brijxs] = fmap mTreeToList [tl, tr, bl, br]
            hshift (i, j, x) = (i + halfs, j, x)
            wshift (i, j, x) = (i, j + halfs, x)
            halfs = s `div` 2

-- utility functions

isZero :: MTree -> Bool
isZero (Zero _) = True
isZero _        = False

ifZeroReplace :: MTree -> MTree
ifZeroReplace tree@(Zero _)   = tree
ifZeroReplace tree@(Leaf _ x) = if x == 0 then Zero 1 else tree
ifZeroReplace tree@(Square s _ _ _ _ _)
                              = if all isZero (subTrees tree)
                                then Zero s else tree

nextPowOf2 :: Integral a => a -> a
nextPowOf2 n = head . dropWhile (< n) $ map (2^) [0..]