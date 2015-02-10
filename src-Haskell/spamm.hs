module SpAMM
( MatrixTree, MatrixList, IndexedMatrixList
, treeTranspose, treeAdd, treeMult, treeMultTol
, combineZeros
, isValidList, isValidIndexedList
, listToIndexedList, indexedListToList
, listToTree, treeToList
, indexedListToTree, treeToIndexedList
, readMatrixList, writeMatrixList
, readIndexedList, writeIndexedList
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

-- other data types

type MatrixList = [[Value]]

type IndexedMatrixList = [((Int, Int), Value)]

isValidList :: MatrixList -> Bool
isValidList list = (not $ null list) && (not $ null (head list)) && (all sameLength list)  
                   where sameLength row = length row == length (head list)

isValidIndexedList :: IndexedMatrixList -> Bool
isValidIndexedList indexedList = (not $ null indexedList) && list == nub list
                                 where list = map fst indexedList

listToIndexedList :: MatrixList -> IndexedMatrixList
listToIndexedList list = filter ((/= 0) . snd) $ zip indexList valueList
                         where rows = length list
                               cols = length $ head list
                               indexList = [(i, j) | i <- [1..rows], j <- [1..cols]]
                               valueList = concat list

indexedListToList :: IndexedMatrixList -> MatrixList
indexedListToList ijxs = map (map value) indexList 
                         where hashTable = Map.fromList ijxs
                               indices = map fst ijxs
                               rows = maximum (map fst indices)
                               cols = maximum (map snd indices)
                               indexList = [[(i, j) | j <- [1..cols]] | i <- [1..rows]]
                               value pair | isNothing check = 0
                                          | otherwise       = fromJust check
                                          where check = Map.lookup pair hashTable

listToTree :: MatrixList -> MatrixTree
listToTree list | isValidList list = snd . fillFromList . addListSize $ list
                | otherwise        = error "Matrix list is invalid"

treeToList :: MatrixTree -> MatrixList
treeToList (Zero t b l r)            = replicate (b - t + 1) $ replicate (r - l + 1) 0
treeToList (Value _ _ _ x)           = [[x]]
treeToList (Row _ _ _ _ ltree rtree) = [concat (treeToList ltree ++ treeToList rtree)]
treeToList (Col _ _ _ _ ttree btree) = treeToList ttree ++ treeToList btree
treeToList (Rect _ _ _ _ _ tltree trtree bltree brtree) =
           zipWith (++) (treeToList tltree ++ treeToList bltree)
                        (treeToList trtree ++ treeToList brtree)

indexedListToTree :: IndexedMatrixList -> MatrixTree
indexedListToTree list | isValidIndexedList list = listToTree . indexedListToList $ list
                       | otherwise               = error "Indexed list is invalid"

treeToIndexedList :: MatrixTree -> IndexedMatrixList
treeToIndexedList = listToIndexedList . treeToList

{-
treeToIndexedList (Zero _ _ _ _)  = []
treeToIndexedList (Value i j _ x) = if x == 0 then [] else [((i,j),x)]
treeToIndexedList (Row _ _ _ _ ltree rtree) = treeToIndexedList ltree ++ 
                                              treeToIndexedList rtree
treeToIndexedList (Col _ _ _ _ ttree btree) = treeToIndexedList ttree ++ 
                                              treeToIndexedList btree
treeToIndexedList (Rect _ _ _ _ _ tltree trtree bltree brtree) =
           treeToIndexedList tltree ++ treeToIndexedList bltree ++
           treeToIndexedList trtree ++ treeToIndexedList brtree
-}

-- reads a MatrixList from a file that lists rows in order on separate lines
-- with entries on each row separated by spaces
readMatrixList :: FilePath -> IO MatrixList
readMatrixList filePath = do contents <- readFile filePath
                             return $ map (map read . words) (lines contents)

-- writes a MatrixList to a file of the format read by readMatrixList
writeMatrixList :: MatrixList -> FilePath -> IO ()
writeMatrixList rowList filePath =
                do handle <- openFile filePath WriteMode
                   mapM_ (hPutStrLn handle) $ map (listToString . (map show)) rowList
                   hClose handle

-- reads an IndexeMatrixList from a file that lists each entry on a separate line
-- with the row index, column index, and value separated by spaces
readIndexedList :: FilePath -> IO IndexedMatrixList
readIndexedList filePath = do contents <- readFile filePath
                              let entryList = map words (lines contents)
                              let splitList = map (splitAt 2) entryList
                              let indices = map (map read . fst) splitList
                              let values = map (read . head . snd) splitList
                              return $ zipWith (\[i,j] x -> ((i,j),x)) indices values

-- writes an IndexedMatrixList to a file of the format read by readIndexedList
writeIndexedList :: IndexedMatrixList -> FilePath -> IO ()
writeIndexedList indexedList filePath =
                 do handle <- openFile filePath WriteMode
                    mapM_ (hPutStrLn handle) $ map pairToString indexedList

-- tree utility functions

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

-- other utility functions

addListSize :: MatrixList -> (MatrixList, Int, Int, Int, Int)
addListSize xss = (xss, 1, length xss, 1, length . head $ xss)

fillFromList :: (MatrixList, Int, Int, Int, Int) -> (Norm, MatrixTree)
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