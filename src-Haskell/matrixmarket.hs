module MatrixMarket
( MatrixTree
, readTreeFromMatrixMarket
, writeTreeToMatrixMarket
) where

-- reads/writes a MatrixTree to MatrixMarket formats

import Data.Char (toLower)
import Data.List (intersperse, nub, transpose)
import qualified Data.Map as Map (fromList, lookup, Map)
import Data.Maybe (fromJust, isNothing)
import MatrixTree (MatrixTree, RowList, rowListToTree, treeToRowList, Value)
import System.IO (hClose, hPutStrLn, openFile, IOMode(WriteMode))

readTreeFromMatrixMarket :: FilePath -> IO MatrixTree
readTreeFromMatrixMarket filePath = readRowListFromMatrixMarket filePath >>=
                                    (return . rowListToTree)

writeTreeToMatrixMarket :: MatrixTree -> String -> FilePath -> IO ()
writeTreeToMatrixMarket tree format filePath =
                        writeRowListToMatrixMarket (treeToRowList tree) format filePath

type ValueList = [Value]

type IndexedList = ((Int,Int),[((Int, Int), Value)])

readRowListFromMatrixMarket :: FilePath -> IO RowList
readRowListFromMatrixMarket filePath = do
           contents <- readFile filePath
           let fileLines = lines contents
           let (first, rest) = (words . fmap toLower $ head fileLines,
                                filter ((/='%') . head) $ tail fileLines)
           if length first < 3 then error $ filePath ++ " header is invalid"
           else case first !! 2 of
                     "array"      -> return $ arrayToRowList rest filePath
                     "coordinate" -> return $ coordinateToRowList rest filePath
                     _            -> error $ filePath ++ " has unrecognized format " ++ (first !! 2)

writeRowListToMatrixMarket :: RowList -> String -> FilePath -> IO ()
writeRowListToMatrixMarket rowList format filePath = do
            let list = case format of "array"      -> rowListToArray rowList
                                      "coordinate" -> rowListToCoordinate rowList
                                      _            -> error $ "Unrecognized format " ++ format
            handle <- openFile filePath WriteMode
            hPutStrLn handle $ "%%MatrixMarket matrix " ++ format ++ " real general"
            mapM_ (hPutStrLn handle) list
            hClose handle

arrayToRowList :: [String] -> FilePath -> RowList
arrayToRowList list filePath | isValid   = transpose $ arrangeInCols values
                             | otherwise = error $ "Array format file " ++ filePath ++
                                                   " is invalid"
                             where isValid = isValidArray list nrows ncols
                                   arrangeInCols [] = []
                                   arrangeInCols list = (take nrows list):
                                                        (arrangeInCols . drop nrows $ list)
                                   [nrows, ncols] = (map read) . words $ head list :: [Int]
                                   values = fmap read $ tail list :: ValueList

coordinateToRowList :: [String] -> FilePath -> RowList
coordinateToRowList list filePath
                    | isValid   = if nonzeros == 0 then replicate (fst size) $ replicate (snd size) 0
                                  else indexedListToRowList indexedList
                    | otherwise = error $ "Coordinate format file " ++ filePath ++ " is invalid"
                    where isValid = isValidCoordinate list firstLine values nonzeros indices size
                          firstLine = fmap read . words $ head list :: [Int]
                          (size, nonzeros) = ((firstLine !! 0, firstLine !! 1), firstLine !! 2)
                          entries = fmap ((splitAt 2) . words) $ tail list
                          indices = fmap (fmap read . fst) entries :: [[Int]]
                          values = fmap (read . head . snd) entries :: [Value]
                          indexedList = (size, zipWith (\[i,j] x -> ((i,j),x)) indices values)

rowListToArray :: RowList -> [String]
rowListToArray rowList = (combineStrings . fmap show $ [length rowList, length $ head rowList]):
                         (fmap show . concat . transpose $ rowList)

rowListToCoordinate :: RowList -> [String]
rowListToCoordinate = indexedListToCoordinate . rowListToIndexedList

indexedListToRowList :: IndexedList -> RowList
indexedListToRowList ((m, n), ijxs)
                     = fmap (fmap value) indices
                       where hashTable = Map.fromList ijxs
                             ijs = fmap fst ijxs
                             indices = [[(i, j) | j <- [1..n]] | i <- [1..m]]
                             value pair | isNothing check = 0
                                        | otherwise       = fromJust check
                                        where check = Map.lookup pair hashTable

rowListToIndexedList :: RowList -> IndexedList
rowListToIndexedList rows = ((nrows, ncols), ijxs)
                            where nrows = length rows
                                  ncols = length $ head rows
                                  indices = [(i, j) | i <- [1..nrows], j <- [1..ncols]]
                                  values = concat rows
                                  ijxs = filter ((/= 0) . snd) $ zip indices values

indexedListToCoordinate :: IndexedList -> [String]
indexedListToCoordinate ((m, n), ijxs) = (combineStrings [show m, show n, show $ length ijxs]):
                                         (fmap pairToString ijxs)

isValidArray :: [String] -> Int -> Int -> Bool
isValidArray list nrows ncols = not (null list) && length (words $ head list) == 2 &&
                                length (tail list) == nrows * ncols

isValidCoordinate :: [String] -> [a] -> [Value] -> Int -> [[Int]] -> (Int,Int) -> Bool
isValidCoordinate list firstLine values nonzeros indices size =
                  not (null list) && length firstLine >= 3 &&
                  all (>= 3) (map (length . words) $ tail list) &&
                  length values == nonzeros && indices == nub indices &&
                  maxCol <= fst size && maxRow <= snd size
                  where [maxCol, maxRow] = if nonzeros == 0 then [1,1] else
                                           fmap maximum [fmap (!! 0) indices, fmap (!! 1) indices]


combineStrings :: [String] -> String
combineStrings = concat . (intersperse " ")

pairToString :: (Show a, Show b, Show c) => ((a,b),c) -> String
pairToString = combineStrings . (\((i,j),x) -> [show i, show j, show x])