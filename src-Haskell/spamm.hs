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

data MatrixTree = Zero  {trow :: Int, brow :: Int, lcol :: Int, rcol :: Int} |
                  Value {row :: Int, col :: Int, normField :: Norm, value :: Value} |
                  Row   {row :: Int, lcol :: Int, rcol :: Int, normField :: Norm,
                         ltree :: MatrixTree, rtree :: MatrixTree} |
                  Col   {col :: Int, trow :: Int, brow :: Int, normField :: Norm,
                         ttree :: MatrixTree, btree :: MatrixTree} |
                  Rect  {trow :: Int, brow :: Int, lcol :: Int, rcol :: Int,
                         normField :: Norm, tltree :: MatrixTree, trtree :: MatrixTree,
                                            bltree :: MatrixTree, brtree :: MatrixTree}
                  deriving (Eq, Show)

-- reading and writing to file

readTreeFromRowList :: FilePath -> IO MatrixTree
readTreeFromRowList filePath = readRowList filePath >>= (return . rowListToTree)

writeTreeToRowList :: MatrixTree -> FilePath -> IO ()
writeTreeToRowList tree filePath = writeRowList (treeToRowList tree) filePath

readTreeFromIndexedList :: FilePath -> IO MatrixTree
readTreeFromIndexedList filePath = readIndexedList filePath >>= (return . indexedListToTree)

writeTreeToIndexedList :: MatrixTree -> FilePath -> IO ()
writeTreeToIndexedList tree filePath =writeIndexedList (treeToIndexedList tree) filePath

-- accessing and calculating norms

norm :: MatrixTree -> Norm
norm (Zero _ _ _ _) = 0
norm tree = normField tree

normNum :: Value -> Norm
normNum = abs

addSubtreeNorms :: [Norm] -> Norm
addSubtreeNorms = sqrt . sum . fmap (^2)

-- matrix algebra

treeTranspose :: MatrixTree -> MatrixTree
treeTranspose (Zero t b l r)            = Zero l r t b
treeTranspose (Value i j x y)           = Value j i x y
treeTranspose (Row i l r x ltree rtree) = Col i l r x (treeTranspose ltree)
                                                      (treeTranspose rtree)
treeTranspose (Col i t b x ttree btree) = Row i t b x (treeTranspose ttree)
                                                      (treeTranspose btree)
treeTranspose (Rect t b l r x tltree trtree bltree brtree) =
               Rect l r t b x (treeTranspose tltree) (treeTranspose bltree)
                              (treeTranspose trtree) (treeTranspose brtree)

treeAdd :: MatrixTree -> MatrixTree -> MatrixTree

treeAdd (Zero i j k l) tree@(Zero m n p q)
        | [i,j,k,l] == [m,n,p,q] = tree
        | otherwise              = error "zero matrices don't match for addition"

treeAdd tree zeroTree@(Zero _ _ _ _) = treeAdd zeroTree tree

treeAdd (Zero i j k l) tree@(Value m n _ _)
        | [i,j,k,l] == [m,m,n,n] = tree
        | otherwise              = error "values don't match for addition"

treeAdd (Zero i j k l) tree@(Row m n p _ _ _)
        | [i,j,k,l] == [m,m,n,p] = tree
        | otherwise              = error "rows don't match for addition"

treeAdd (Zero i j k l) tree@(Col m n p _ _ _)
        | [i,j,k,l] == [n,p,m,m] = tree
        | otherwise              = error "columns don't match for addition"

treeAdd (Zero i j k l) tree@(Rect m n p q _ _ _ _ _)
        | [i,j,k,l] == [m,n,p,q] = tree
        | otherwise              = error "rectangles don't match for addition"

treeAdd (Value i j _ x) (Value k l _ y)
        | [i,j] == [k,l] = if x + y == 0 then Zero i i j j
                           else Value i j (normNum $ x + y) (x + y)
        | otherwise      = error "values don't match for addition"

treeAdd (Row i l r _ ltree1 rtree1) (Row j m s _ ltree2 rtree2)
        | [i,l,r] == [j,m,s] = ifZeroReplace (Row i l r x lsum rsum)
        | otherwise          = error "rows don't match for addition"
        where lsum = ltree1 `treeAdd` ltree2
              rsum = rtree1 `treeAdd` rtree2
              x    = addSubtreeNorms [norm lsum, norm rsum]

treeAdd (Col i t b _ ttree1 btree1) (Col j u c _ ttree2 btree2)
        | [i,t,b] == [j,u,c] = ifZeroReplace (Col i t b x tsum bsum)
        | otherwise          = error "columns don't match for addition"
        where tsum = ttree1 `treeAdd` ttree2
              bsum = btree1 `treeAdd` btree2
              x    = addSubtreeNorms [norm tsum, norm bsum]

treeAdd (Rect t b l r _ tltree1 trtree1 bltree1 brtree1)
        (Rect u c m s _ tltree2 trtree2 bltree2 brtree2)
        | [t,b,l,r] == [u,c,m,s] = ifZeroReplace (Rect t b l r x tlsum trsum blsum brsum)
        | otherwise              = error "rectangles don't match for addition"
        where tlsum = tltree1 `treeAdd` tltree2
              trsum = trtree1 `treeAdd` trtree2
              blsum = bltree1 `treeAdd` bltree2
              brsum = brtree1 `treeAdd` brtree2
              x     = addSubtreeNorms . fmap norm $ [tlsum, trsum, blsum, brsum]

treeAdd _ _ = error "shapes don't match for addition"

treeMult :: MatrixTree -> MatrixTree -> MatrixTree
tree1 `treeMult` tree2 = treeMultTol tree1 tree2 0

treeMultTol :: MatrixTree -> MatrixTree -> Double -> MatrixTree

treeMultTol (Zero i j k l) (Zero m n p q) _
            | [k,l] == [m,n] = Zero i j p q
            | otherwise      = error "zero matrices don't match for multiplication"

treeMultTol tree zeroTree@(Zero _ _ _ _) _ =
            treeTranspose $ treeMultTol (treeTranspose zeroTree) (treeTranspose tree) 0

treeMultTol (Zero i j k l) (Value m n _ _) _
            | [k,l] == [m,m] = Zero i j n n
            | otherwise      = error "values don't match for multiplication"

treeMultTol (Zero i j k l) (Row m n p _ _ _) _
            | [k,l] == [m,m] = Zero i j n p
            | otherwise      = error "rows don't match for multiplication"

treeMultTol (Zero i j k l) (Col m n p _ _ _) _
            | [k,l] == [n,p] = Zero i j m m
            | otherwise      = error "columns don't match for multiplication"

treeMultTol (Zero i j k l) (Rect m n p q _ _ _ _ _) _
            | [k,l] == [m,n] = Zero i j p q
            | otherwise      = error "rectangles don't match for multiplication"

--only this treeMultTol Col Value or the one below can be uncommented
treeMultTol tree@(Col _ _ _ _ _ _) val@(Value _ _ _ _) tol =
            treeTranspose $ treeMultTol (treeTranspose val) (treeTranspose tree) tol

-- only this treeMultTol Rect Col or the one below can be uncommented
treeMultTol tree1@(Rect _ _ _ _ _ _ _ _ _) tree2@(Col _ _ _ _ _ _) tol =
            treeTranspose $ treeMultTol (treeTranspose tree2) (treeTranspose tree1) tol

treeMultTol (Value i j m x) (Value k l n y) tol
            | j == k    = if m * n <= tol then Zero i i l l
                          else Value i l (normNum $ x * y) (x * y)
            | otherwise = error "values don't match for multiplication"

treeMultTol val@(Value i j m x) (Row k l r n ltree rtree) tol
            | j == k    = if m * n <= tol then Zero i i l r
                          else ifZeroReplace (Row i l r x lmult rmult)
            | otherwise = error "value and row don't match for multiplication"
            where lmult = treeMultTol val ltree tol
                  rmult = treeMultTol val rtree tol
                  x     = addSubtreeNorms [norm lmult, norm rmult]

treeMultTol (Row i l r m ltree rtree) (Col j t b n ttree btree) tol
            | [l,r] == [t,b] = if m * n <= tol then Zero i i j j
                               else treeAdd (treeMultTol ltree ttree tol)
                                            (treeMultTol rtree btree tol)
            | otherwise      = error "row and column don't match for multiplication"

treeMultTol (Row k l1 r1 m ltree rtree) (Rect t b l2 r2 n tltree trtree bltree brtree) tol
            | [l1,r1] == [t,b] = if m * n <= tol then Zero k k l2 r2
                                 else ifZeroReplace (Row k l2 r2 x lmult rmult)
            | otherwise        = error "row and rectangle don't match for multiplication"
            where lmult = treeAdd (treeMultTol ltree tltree tol)
                                  (treeMultTol rtree bltree tol)
                  rmult = treeAdd (treeMultTol ltree trtree tol)
                                  (treeMultTol rtree brtree tol)
                  x     = addSubtreeNorms [norm lmult, norm rmult]

-- only this treeMultTol Col Value or the one above can be uncommented
{-
treeMultTol (Col i t b m ttree btree) val@(Value j k n x) tol
            | i == j    = if m * n <= tol then Zero t b i i
                          else ifZeroReplace (Col i t b x tmult bmult)
            | otherwise = error "column and value don't match for multiplication"
            where tmult = treeMultTol ttree val tol
                  bmult = treeMultTol btree val tol
                  x     = addSubtreeNorms [norm tmult, norm bmult]
-}

treeMultTol (Col i t b m ttree btree) (Row k l r n ltree rtree) tol
            | [t,b] == [l,r] = if m * n <= tol then Zero t b l r
                               else ifZeroReplace (Rect t b l r x
                                                        tlmult trmult blmult brmult)
            | otherwise      = error "column and row don't match for multiplication"
            where tlmult = treeMultTol ttree ltree tol
                  trmult = treeMultTol ttree rtree tol
                  blmult = treeMultTol btree ltree tol
                  brmult = treeMultTol btree rtree tol
                  x      = addSubtreeNorms . fmap norm $ [tlmult, trmult, blmult, brmult]

-- only this treeMultTol Rect Col or the one above can be uncommented
{-
treeMultTol (Rect t1 b1 l r m tltree trtree bltree brtree) (Col i t2 b2 n ttree btree) tol
            | [l,r] == [t2,b2] = if m * n <= tol then Zero t1 b1 i i
                                 else ifZeroReplace (Col i t1 b1 x tmult bmult)
            | otherwise        = error "rectangle and column don't match for multiplication"
            where tmult = treeAdd (treeMultTol tltree ttree tol)
                                  (treeMultTol trtree btree tol)
                  bmult = treeAdd (treeMultTol bltree ttree tol)
                                  (treeMultTol brtree btree tol)
                  x     = addSubtreeNorms [norm tmult, norm bmult]
-}

treeMultTol (Rect t1 b1 l1 r1 m tltree1 trtree1 bltree1 brtree1)
            (Rect t2 b2 l2 r2 n tltree2 trtree2 bltree2 brtree2) tol
            | [l1,r1] == [t2,b2] = if m * n <= tol then Zero t1 b1 l2 r2
                                   else ifZeroReplace (Rect t1 b1 l2 r2 x
                                                            tlmult trmult blmult brmult)
            | otherwise          = error "rectangles don't match for multiplication"
            where tlmult = treeAdd (treeMultTol tltree1 tltree2 tol)
                                   (treeMultTol trtree1 bltree2 tol)
                  trmult = treeAdd (treeMultTol tltree1 trtree2 tol)
                                   (treeMultTol trtree1 brtree2 tol)
                  blmult = treeAdd (treeMultTol bltree1 tltree2 tol)
                                   (treeMultTol brtree1 bltree2 tol)
                  brmult = treeAdd (treeMultTol bltree1 trtree2 tol)
                                   (treeMultTol brtree1 brtree2 tol)
                  x      = addSubtreeNorms . fmap norm $ [tlmult, trmult, blmult, brmult]

treeMultTol _ _ _ = error "shapes don't match for multiplication"

-- utility functions

isZero :: MatrixTree -> Bool
isZero (Zero _ _ _ _) = True
isZero _              = False

ifZeroReplace :: MatrixTree -> MatrixTree
ifZeroReplace tree@(Zero _ _ _ _) = tree
ifZeroReplace tree@(Value i j _ x) = if x == 0 then Zero i i j j else tree
ifZeroReplace tree@(Row i l r _ ltree rtree) =
              if isZero ltree && isZero rtree then Zero i i l r else tree
ifZeroReplace tree@(Col i t b _ ttree btree) =
              if isZero ttree && isZero btree then Zero t b i i else tree
ifZeroReplace tree@(Rect t b l r _ tltree trtree bltree brtree) =
              if isZero tltree && isZero trtree && isZero bltree && isZero brtree
              then Zero t b l r else tree

combineZeros :: MatrixTree -> MatrixTree
combineZeros (Row i l r _ ltree rtree) =
             ifZeroReplace (Row i l r x newltree newrtree)
             where newltree = combineZeros ltree
                   newrtree = combineZeros rtree
                   x = addSubtreeNorms [norm newltree, norm newrtree]
combineZeros (Col i t b _ ttree btree) =
             ifZeroReplace (Col i t b x newttree newbtree)
             where newttree = combineZeros ttree
                   newbtree = combineZeros btree
                   x = addSubtreeNorms [norm newttree, norm newbtree]
combineZeros (Rect t b l r _ tltree trtree bltree brtree) =
             ifZeroReplace (Rect t b l r x newtltree newtrtree newbltree newbrtree)
             where newtltree = combineZeros tltree
                   newtrtree = combineZeros trtree
                   newbltree = combineZeros bltree
                   newbrtree = combineZeros brtree
                   x = addSubtreeNorms . fmap norm $ [newtltree, newtrtree,
                                                      newbltree, newbrtree]
combineZeros tree = ifZeroReplace tree

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
                                       ijs = map fst ijxs
                                       maxRow = maximum $ map fst ijs
                                       maxCol = maximum $ map snd ijs
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
                     = map (map value) indices
                       where hashTable = Map.fromList ijxs
                             ijs = map fst ijxs
                             indices = [[(i, j) | j <- [1..n]] | i <- [1..m]]
                             value pair | isNothing check = 0
                                        | otherwise       = fromJust check
                                        where check = Map.lookup pair hashTable

rowListToTree :: RowList -> MatrixTree
rowListToTree rows | isValidRowList rows = snd . fillFromList . addListSize $ rows
                   | otherwise           = error "Row list is invalid"

treeToRowList :: MatrixTree -> RowList
treeToRowList (Zero t b l r)            = replicate (b - t + 1) $ replicate (r - l + 1) 0
treeToRowList (Value _ _ _ x)           = [[x]]
treeToRowList (Row _ _ _ _ ltree rtree) = [concat (treeToRowList ltree ++ treeToRowList rtree)]
treeToRowList (Col _ _ _ _ ttree btree) = treeToRowList ttree ++ treeToRowList btree
treeToRowList (Rect _ _ _ _ _ tltree trtree bltree brtree) =
              zipWith (++) (treeToRowList tltree ++ treeToRowList bltree)
                           (treeToRowList trtree ++ treeToRowList brtree)

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
                          return $ map (map read . words) (lines contents)

-- writes a RowList to a file of the format read by readRowList
writeRowList :: RowList -> FilePath -> IO ()
writeRowList rowList filePath =
             do handle <- openFile filePath WriteMode
                mapM_ (hPutStrLn handle) $ map (listToString . (map show)) rowList
                hClose handle

-- reads an IndexeMatrixList from a file that lists each entry on a separate line
-- with the row index, column index, and value separated by spaces
readIndexedList :: FilePath -> IO IndexedList
readIndexedList filePath =
                do contents <- readFile filePath
                   let filteredContents = map words . filter ((/= '%') . head) $
                                          lines contents
                   let size = (\[x,y] -> (x,y)) . (map read) . take 2 $ head filteredContents
                   let entries = map (splitAt 2) $ tail filteredContents
                   let indices = map (map read . fst) entries
                   let values = map (read . head . snd) entries
                   return (size, zipWith (\[i,j] x -> ((i,j),x)) indices values)

-- writes an IndexedList to a file of the format read by readIndexedList
writeIndexedList :: IndexedList -> FilePath -> IO ()
writeIndexedList indexedList filePath =
                 do let (size, ijxs) = indexedList
                    handle <- openFile filePath WriteMode
                    hPutStrLn handle $ (listToString . map show) [fst size, snd size]
                    mapM_ (hPutStrLn handle) $ map pairToString ijxs

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
addListSize xss = (xss, 1, length xss, 1, length . head $ xss)

fillFromList :: (RowList, Int, Int, Int, Int) -> (Norm, MatrixTree)
fillFromList ([[x]], i, _, j, _) = if x == 0 then (0, Zero i i j j)
                                   else (normNum x, Value i j (normNum x) x)
fillFromList (xss, t, b, l, r)
             | all0 xss  = (0, Zero t b l r)
             | t == b    = (rowNorm, Row t l r rowNorm (tree lrow) (tree rrow))
             | l == r    = (colNorm, Col l t b colNorm (tree tcol) (tree bcol))
             | otherwise = (rectNorm, Rect t b l r rectNorm (tree re11) (tree re12)
                                                            (tree re21) (tree re22))
             where all0 = and . fmap (== 0) . concat
                   rowNorm = calcNorm [lrow, rrow]
                   colNorm = calcNorm [tcol, bcol]
                   rectNorm = calcNorm [re11, re12, re21, re22]
                   calcNorm = addSubtreeNorms . fmap (fst . fillFromList)
                   lrow = (lxss, t, t, l, midpt l r)
                   rrow = (rxss, t, t, midpt l r + 1, r)
                   tcol = (txss, t, midpt t b, l, l)
                   bcol = (bxss, midpt t b + 1, b, l, l)
                   re11 = (xss11, t, midpt t b, l, midpt l r)
                   re12 = (xss12, t, midpt t b, midpt l r + 1, r)
                   re21 = (xss21, midpt t b + 1, b, l, midpt l r)
                   re22 = (xss22, midpt t b + 1, b, midpt l r + 1, r)
                   (lxss, rxss) = halveRowList xss
                   (txss, bxss) = halveColList xss
                   [xss11, xss12, xss21, xss22] = partitionList xss
                   midpt m n = m + (n - m - 1) `div` 2
                   tree = snd . fillFromList

halveList :: [a] -> ([a],[a])
halveList xs = splitAt (length xs `div` 2) xs

halveRowList :: [[a]] -> ([[a]],[[a]])
halveRowList xss = let (lxs, rxs) = head $ fmap halveList xss in ([lxs], [rxs])

halveColList :: [[a]] -> ([[a]],[[a]])
halveColList xss = (transpose txss, transpose bxss)
                   where (txss, bxss) = halveRowList . transpose $ xss

partitionList :: [[a]] -> [[[a]]]
partitionList xss = fmap transpose [tlxss, trxss, blxss, brxss]
                    where (txss, bxss)   = halveList xss
                          (tlxss, trxss) = halveList . transpose $ txss
                          (blxss, brxss) = halveList . transpose $ bxss

listToString :: [String] -> String
listToString = concat . (intersperse " ")

pairToString :: (Show a, Show b, Show c) => ((a,b),c) -> String
pairToString = listToString . (\((i,j),x) -> [show i, show j, show x])