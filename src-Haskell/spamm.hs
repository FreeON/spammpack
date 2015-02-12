module SpAMM
( MatrixTree
, treeTranspose, treeAdd, treeMult, treeMultTol
, combineZeros
, readTreeFromRowList, writeTreeToRowList
, readTreeFromIndexedList, writeTreeToIndexedList
) where

import Data.List (intersperse, nub, transpose)
import qualified Data.Map as Map (fromList, lookup, Map)
import Data.Maybe (fromJust, isNothing)
import System.IO (hClose, hPutStrLn, openFile, IOMode(WriteMode))

-- data structure

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

-- reading from and writing to a file

readTreeFromRowList :: FilePath -> IO MatrixTree
readTreeFromRowList filePath = readRowList filePath >>= (return . rowListToTree)

writeTreeToRowList :: MatrixTree -> FilePath -> IO ()
writeTreeToRowList tree filePath = writeRowList (treeToRowList tree) filePath

readTreeFromIndexedList :: FilePath -> IO MatrixTree
readTreeFromIndexedList filePath = readIndexedList filePath >>= (return . indexedListToTree)

writeTreeToIndexedList :: MatrixTree -> FilePath -> IO ()
writeTreeToIndexedList tree filePath =writeIndexedList (treeToIndexedList tree) filePath

-- accessing and calculating norms

getNorm :: MatrixTree -> Norm
getNorm (Zero _ _ _ _) = 0
getNorm tree = norm tree

valueNorm :: Value -> Norm
valueNorm = abs

addSubtreeNorms :: [Norm] -> Norm
addSubtreeNorms = sqrt . sum . fmap (^2)

-- matrix algebra

treeTranspose :: MatrixTree -> MatrixTree
treeTranspose (Zero t l h w)               = Zero l t w h
treeTranspose (Value i j x y)              = Value j i x y
treeTranspose (Rect t l h w x tl tr bl br) = Rect l t w h x
                                                  (treeTranspose tl) (treeTranspose bl)
                                                  (treeTranspose tr) (treeTranspose br)

treeAdd :: MatrixTree -> MatrixTree -> MatrixTree

treeAdd (Zero i j k l) tree@(Zero m n p q)
        | [i,j,k,l] == [m,n,p,q] = tree
        | otherwise              = error "zero matrices don't match for addition"

treeAdd tree zeroTree@(Zero _ _ _ _) = treeAdd zeroTree tree

treeAdd (Zero t l h w) tree@(Value i j _ _)
        | [t,l,h,w] == [i,j,1,1] = tree
        | otherwise              = error "values don't match for addition"

treeAdd (Zero t l h w) tree@(Rect u m j v _ _ _ _ _)
        | [t,l,h,w] == [u,m,j,v] = tree
        | otherwise              = error "rectangles don't match for addition"

treeAdd (Value i j _ x) (Value k l _ y)
        | [i,j] == [k,l] = if x + y == 0 then Zero i j 1 1
                           else Value i j (valueNorm $ x + y) (x + y)
        | otherwise      = error "values don't match for addition"

treeAdd rect1@(Rect t l h w _ _ _ _ _) rect2@(Rect u m j v _ _ _ _ _)
        | [t,l,h,w] == [u,m,j,v] = ifZeroReplace (Rect t l h w x tlsum trsum blsum brsum)
        | otherwise              = error "rectangles don't match for addition"
        where Rect _ _ _ _ _ tl1 tr1 bl1 br1 = rectOrder rect1
              Rect _ _ _ _ _ tl2 tr2 bl2 br2 = rectOrder rect2
              tlsum = tl1 `treeAdd` tl2 ; trsum = tr1 `treeAdd` tr2
              blsum = bl1 `treeAdd` bl2 ; brsum = br1 `treeAdd` br2
              x = addSubtreeNorms . fmap getNorm $ [tlsum, trsum, blsum, brsum]

treeAdd _ _ = error "shapes don't match for addition"

treeMult :: MatrixTree -> MatrixTree -> MatrixTree
tree1 `treeMult` tree2 = treeMultTol tree1 tree2 0

treeMultTol :: MatrixTree -> MatrixTree -> Double -> MatrixTree

treeMultTol (Zero t1 l1 h1 w1) (Zero t2 l2 h2 w2) _
            | [l1,w1] == [t2,h2] = Zero t1 l2 h1 w2
            | otherwise          = error "zero matrices don't match for multiplication"

treeMultTol tree zeroTree@(Zero _ _ _ _) _ =
            treeTranspose $ treeMult (treeTranspose zeroTree) (treeTranspose tree)

treeMultTol (Zero t l h w) (Value i j _ _) _
            | l == i && w `elem` [0,1] = Zero t j h w
            | otherwise      = error "zero and value don't match for multiplication"

treeMultTol (Zero t1 l1 h1 w1) (Rect t2 l2 h2 w2 _ _ _ _ _) _
            | [l1,w1] == [t2,h2] = Zero t1 l2 h1 w2
            | otherwise          = error "rectangles don't match for multiplication"

treeMultTol tree1@(Rect _ _ _ _ _ _ _ _ _) tree2@(Value _ _ _ _) tol =
            treeTranspose $ treeMultTol (treeTranspose tree2) (treeTranspose tree1) tol

treeMultTol (Value i j m x) (Value k l n y) tol
            | j == k    = if m * n <= tol then Zero i l 1 1
                          else Value i l (valueNorm $ x * y) (x * y)
            | otherwise = error "values don't match for multiplication"

treeMultTol val@(Value i j m _) (Rect t l h w n tl tr bl br) tol
            | [j,1] == [t,h] = if m * n <= tol then Zero i l 1 w
                               else ifZeroReplace (Rect i l 1 w x
                                                        tlmult trmult blmult brmult)
            | otherwise      = error "value and rectangle don't match for multiplication"
            where [tlmult, trmult, blmult, brmult] = fmap multiply [tl, tr, bl, br]
                  multiply tree = treeMultTol val tree tol
                  x = addSubtreeNorms . fmap getNorm $ [tlmult, trmult, blmult, brmult]

treeMultTol (Rect t1 l1 h1 w1 m tl1 tr1 bl1 br1) (Rect t2 l2 h2 w2 n tl2 tr2 bl2 br2) tol
            | [l1,w1] == [t2,h2] = if m * n <= tol then Zero t1 l2 h1 w2
                                   else ifZeroReplace (Rect t1 l2 h1 w2 x
                                                            tlmult trmult blmult brmult)
            | otherwise          = error "rectangles don't match for multiplication"
            where tlmult = treeAdd (treeMultTol tl1 tl2 tol) (treeMultTol tr1 bl2 tol)
                  trmult = treeAdd (treeMultTol tl1 tr2 tol) (treeMultTol tr1 br2 tol)
                  blmult = treeAdd (treeMultTol bl1 tl2 tol) (treeMultTol br1 bl2 tol)
                  brmult = treeAdd (treeMultTol bl1 tr2 tol) (treeMultTol br1 br2 tol)
                  x = addSubtreeNorms . fmap getNorm $ [tlmult, trmult, blmult, brmult]

-- utility functions

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

-- other data types

type RowList = [[Value]]

type IndexedList = ((Int,Int),[((Int, Int), Value)])

isValidRowList :: RowList -> Bool
isValidRowList rows = not (null rows) && not (null $ head rows) && all sameLength rows
                      where sameLength row = length row == length (head rows)

isValidIndexedList :: IndexedList -> Bool
isValidIndexedList indexedList = not (null ijxs)    && ijs == nub ijs &&
                                 maxRow <= fst size && maxCol <= snd size
                                 where ijxs = snd indexedList
                                       ijs = fmap fst ijxs
                                       maxRow = maximum $ fmap fst ijs
                                       maxCol = maximum $ fmap snd ijs
                                       size = fst indexedList

rowListToIndexedList :: RowList -> IndexedList
rowListToIndexedList rows = ((nrows, ncols), ijxs)
                            where nrows = length rows
                                  ncols = length $ head rows
                                  indices = [(i, j) | i <- [1..nrows], j <- [1..ncols]]
                                  values = concat rows
                                  ijxs = filter ((/= 0) . snd) $ zip indices values

indexedListToRowList :: IndexedList -> RowList
indexedListToRowList ((m, n), ijxs)
                     = fmap (fmap value) indices
                       where hashTable = Map.fromList ijxs
                             ijs = fmap fst ijxs
                             indices = [[(i, j) | j <- [1..n]] | i <- [1..m]]
                             value pair | isNothing check = 0
                                        | otherwise       = fromJust check
                                        where check = Map.lookup pair hashTable

rowListToTree :: RowList -> MatrixTree
rowListToTree rows | isValidRowList rows = snd . fillFromList . addListSize $ rows
                   | otherwise           = error "Row list is invalid"

treeToRowList :: MatrixTree -> RowList
treeToRowList (Zero t l h w)               = replicate h $ replicate w 0
treeToRowList (Value _ _ _ x)              = [[x]]
treeToRowList (Rect _ _ _ _ _ tl tr bl br) =
              zipWith (++) (treeToRowList tl ++ treeToRowList bl)
                           (treeToRowList tr ++ treeToRowList br)

indexedListToTree :: IndexedList -> MatrixTree
indexedListToTree list
                  | isValidIndexedList list = rowListToTree . indexedListToRowList $ list
                  | otherwise               = error "Indexed list is invalid"

treeToIndexedList :: MatrixTree -> IndexedList
treeToIndexedList = rowListToIndexedList . treeToRowList

-- reads a RowList from a file that lists rows in order on separate lines
-- with entries on each row separated by spaces
readRowList :: FilePath -> IO RowList
readRowList filePath = do contents <- readFile filePath
                          return $ fmap (fmap read . words) (lines contents)

-- writes a RowList to a file of the format read by readRowList
writeRowList :: RowList -> FilePath -> IO ()
writeRowList rowList filePath =
             do handle <- openFile filePath WriteMode
                mapM_ (hPutStrLn handle) $ fmap (listToString . (fmap show)) rowList
                hClose handle

-- reads an IndexeMatrixList from a file that lists each entry on a separate line
-- with the row index, column index, and value separated by spaces
readIndexedList :: FilePath -> IO IndexedList
readIndexedList filePath =
                do contents <- readFile filePath
                   let filteredContents = fmap words . filter ((/= '%') . head) $
                                          lines contents
                   let size = (\[x,y] -> (x,y)) . (fmap read) $ head filteredContents
                   let entries = fmap (splitAt 2) $ tail filteredContents
                   let indices = fmap (fmap read . fst) entries
                   let values = fmap (read . head . snd) entries
                   return (size, zipWith (\[i,j] x -> ((i,j),x)) indices values)

-- writes an IndexedList to a file of the format read by readIndexedList
writeIndexedList :: IndexedList -> FilePath -> IO ()
writeIndexedList indexedList filePath =
                 do let (size, ijxs) = indexedList
                    handle <- openFile filePath WriteMode
                    hPutStrLn handle $ (listToString . fmap show) [fst size, snd size]
                    mapM_ (hPutStrLn handle) $ fmap pairToString ijxs

-- tests of internal functions

testRowIndexedEq :: IO ()
testRowIndexedEq = print $ rowListToTree testRowList == indexedListToTree testIndexedList

testRowList = [[0, 0, 0, 7,  0, 0],
               [0, 0, 0, 0,  0, 0],
               [0, 0, 0, 0,  0, 0],
               [0, 0, 0, 0, 12, 0],
               [0, 0, 0, 0,  0, 0],
               [0, 0, 0, 0,  0, 0],
               [0, 0, 0, 0,  0, 0],
               [0, 0, 0, 0,  0, 0],
               [0, 3, 0, 0,  0, 0],
               [0, 0, 0, 0,  0, 0]] :: RowList

testIndexedList = ( (10, 6), [((1, 4), 7),
                              ((9, 2), 3),
                              ((4, 5), 12)] ) :: IndexedList

-- other utility functions

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

halveList :: [a] -> ([a],[a])
halveList xs = splitAt (length xs `div` 2) xs

partitionList :: [[a]] -> [[[a]]]
partitionList xss = fmap transpose [tlxss, trxss, blxss, brxss]
                    where (txss, bxss)   = halveList xss
                          (tlxss, trxss) = halveList . transpose $ txss
                          (blxss, brxss) = halveList . transpose $ bxss

listToString :: [String] -> String
listToString = concat . (intersperse " ")

pairToString :: (Show a, Show b, Show c) => ((a,b),c) -> String
pairToString = listToString . (\((i,j),x) -> [show i, show j, show x])