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
           if null (words contents) then error $ filePath ++ " is empty"
           else do
            let fileLines = lines contents
            let first = words . head $ fileLines
            if length first < 4 then error $ filePath ++ " header is invalid"
            else do
             let rest = filter (not . null) . fmap words . filter ((/='%') . head) $ tail fileLines
             if null rest then error $ filePath ++ " has no contents"
             else do
              let command =
                   case fmap toLower (first !! 2) of
                        "array"      -> arrayToRowList
                        "coordinate" -> coordinateToRowList
                        _            -> error $ filePath ++ " has unrecognized format " ++ (first !! 2)
              return $ command rest filePath

writeRowListToMatrixMarket :: RowList -> String -> FilePath -> IO ()
writeRowListToMatrixMarket rowList format filePath = do
            let list = case format of "array"      -> rowListToArray rowList
                                      "coordinate" -> rowListToCoordinate rowList
                                      _            -> error $ "Unrecognized format " ++ format
            handle <- openFile filePath WriteMode
            hPutStrLn handle $ "%%MatrixMarket matrix " ++ format ++ " real general"
            mapM_ (hPutStrLn handle) list
            hClose handle

arrayToRowList :: [[String]] -> FilePath -> RowList
arrayToRowList lineList filePath
               | isValid   = transpose $ arrangeInCols values
               | otherwise = error $ "Array format file " ++ filePath ++ " is invalid"
               where isValid = length firstLine == 2 && all ((== 1) . length) valueLines &&
                               all (not . null) size && length valueList == nrows * ncols &&
                               all (not . null) valueReads
                     (firstLine, valueLines) = (head lineList, tail lineList)
                     size = fmap reads firstLine :: [[(Int,String)]]
                     [nrows, ncols] = fmap (fst . head) size
                     valueList = concat valueLines
                     valueReads = fmap reads valueList :: [[(Value,String)]]
                     values = fmap (fst . head) valueReads
                     arrangeInCols [] = []
                     arrangeInCols list = (take nrows list):
                                          (arrangeInCols . drop nrows $ list)

coordinateToRowList :: [[String]] -> FilePath -> RowList
coordinateToRowList lineList filePath
                    | isValid   = if nonzeros == 0 then replicate (fst size) $ replicate (snd size) 0
                                  else indexedListToRowList indexedList
                    | otherwise = error $ "Coordinate format file " ++ filePath ++ " is invalid"
                    where isValid = all ((== 3) . length) lineList && all (not . null) firstLineNums &&
                                    length entryLines == nonzeros &&
                                    all (all (not . null)) indexReads && indices == nub indices &&
                                    maxCol <= fst size && maxRow <= snd size &&
                                    all (not . null) valueReads
                          (firstLine, entryLines) = (head lineList, tail lineList)
                          firstLineNums = fmap reads firstLine :: [[(Int,String)]]
                          firstLineVals = fmap (fst . head) firstLineNums
                          size = (firstLineVals !! 0, firstLineVals !! 1)
                          nonzeros = firstLineVals !! 2
                          splitEntries = fmap (splitAt 2) entryLines
                          indexLines = fmap fst splitEntries
                          indexReads = fmap (fmap reads) indexLines :: [[[(Int,String)]]]
                          indices = fmap (fmap (fst . head)) indexReads
                          [maxCol, maxRow] = if nonzeros == 0 then [1,1] else
                                             fmap maximum [fmap (!! 0) indices, fmap (!! 1) indices]
                          valueLines = concat . fmap snd $ splitEntries
                          valueReads = fmap reads valueLines :: [[(Value,String)]]
                          values = fmap (fst . head) valueReads
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

combineStrings :: [String] -> String
combineStrings = concat . (intersperse " ")

pairToString :: (Show a, Show b, Show c) => ((a,b),c) -> String
pairToString = combineStrings . (\((i,j),x) -> [show i, show j, show x])