module MatrixMarket
( MatrixList
, readFromMatrixMarket
, writeToMatrixMarket
) where

-- reads/writes MatrixMarket files to a generic matrix data type

import Data.Char (toLower)
import Data.List (intersperse, nub)
import qualified Data.Map as Map (fromList, lookup, Map)
import Data.Maybe (fromJust, isNothing)
import System.IO (hClose, hPutStrLn, openFile, IOMode(WriteMode))

type MatrixList = ( Int, Int, [ (Int, Int, Double) ] )

-- matrix height, width, and list of entries (row, column, value) ;
-- including zero-value entries is optional

readFromMatrixMarket :: FilePath -> IO MatrixList
readFromMatrixMarket filePath = readFile filePath >>=
                                (return . matrixMarketToMatrixList filePath)

writeToMatrixMarket :: MatrixList -> String -> FilePath -> IO ()
writeToMatrixMarket matrixList format filePath = do
            let list = case format of "array"      -> matrixListToArray matrixList
                                      "coordinate" -> matrixListToCoordinate matrixList
                                      _            -> error $ "Unrecognized format " ++ format
            handle <- openFile filePath WriteMode
            hPutStrLn handle $ "%%MatrixMarket matrix " ++ format ++ " real general"
            mapM_ (hPutStrLn handle) list
            hClose handle

matrixMarketToMatrixList :: FilePath -> String -> MatrixList
matrixMarketToMatrixList filePath contents
      | null (words contents) = error $ filePath ++ " is empty"
      | length first < 4      = error $ filePath ++ " header is invalid"
      | null rest             = error $ filePath ++ " has no contents"
      | otherwise             = command rest filePath
      where fileLines = lines contents
            first = words . head $ fileLines
            rest = filter ((/='%') . head . head) . filter (not . null) . fmap words $ tail fileLines
            command = case fmap toLower (first !! 2) of
                      "array"      -> arrayToMatrixList
                      "coordinate" -> coordinateToMatrixList
                      _            -> error $ filePath ++ " has unrecognized format " ++ (first !! 2)

matrixListToArray :: MatrixList -> [String]
matrixListToArray (m, n, ijxs) = (combineStrings [show m, show n]):showVals
                   where pairList = fmap (\(i, j, x) -> ((i,j), x)) ijxs
                         hashTable = Map.fromList pairList
                         indices = [(i, j) | j <- [1..n], i <- [1..m]]
                         showVals = fmap (show . value) indices
                         value pair | isNothing check = 0
                                    | otherwise       = fromJust check
                                    where check = Map.lookup pair hashTable

matrixListToCoordinate :: MatrixList -> [String]
matrixListToCoordinate (m, n, ijxs) = (combineStrings [show m, show n, show $ length ijxs]):
                                         (fmap pairToString ijxs)

combineStrings :: [String] -> String
combineStrings = concat . (intersperse " ")

pairToString :: (Show a, Show b, Show c) => (a,b,c) -> String
pairToString = combineStrings . (\(i,j,x) -> [show i, show j, show x])

arrayToMatrixList :: [[String]] -> FilePath -> MatrixList
arrayToMatrixList lineList filePath
     | length firstLine /= 2             = error $ filePath ++ " matrix size invalid"
     | any ((/= 1) . length) valueLines  = error $ filePath ++ " has invalid value lines"
     | any null size                     = error $ filePath ++ " matrix size invalid"
     | length valueList /= nrows * ncols = error $ filePath ++ " has wrong matrix size"
     | any null valueReads               = error $ filePath ++ " has invalid values"
     | otherwise                         = (nrows, ncols, ijxs)
     where (firstLine, valueLines) = (head lineList, tail lineList)
           size = fmap reads firstLine :: [[(Int,String)]]
           [nrows, ncols] = fmap (fst . head) size
           valueList = concat valueLines
           valueReads = fmap reads valueList :: [[(Double,String)]]
           values = fmap (fst . head) valueReads
           indices = [(i, j) | j <- [1..ncols], i <- [1..nrows]]
           ijxs = fmap (\((i, j), x) -> (i, j, x)) . filter ((/= 0) . snd) . zip indices $ values

coordinateToMatrixList :: [[String]] -> FilePath -> MatrixList
coordinateToMatrixList lineList filePath
                    | any ((/= 3) . length) lineList   = error $ filePath ++ " has wrong-length line"
                    | any null firstLineNums           = error $ filePath ++ " has invalid size line"
                    | length entryLines /= nonzeros    = error $ filePath ++ " has wrong no. of nonzeros"
                    | any (any null) indexReads        = error $ filePath ++ " has invalid indices"
                    | indices /= nub indices           = error $ filePath ++ " has duplicate indices"
                    | maxRow > nrows || maxCol > ncols = error $ filePath ++ " has indices outside range"
                    | any null valueReads              = error $ filePath ++ " has invalid values"
                    | otherwise                        = (nrows, ncols, if nonzeros == 0 then [] else ijxs)
                    where (firstLine, entryLines) = (head lineList, tail lineList)
                          firstLineNums = fmap reads firstLine :: [[(Int,String)]]
                          firstLineVals = fmap (fst . head) firstLineNums
                          [nrows, ncols] = take 2 firstLineVals
                          nonzeros = firstLineVals !! 2
                          splitEntries = fmap (splitAt 2) entryLines
                          indexLines = fmap fst splitEntries
                          indexReads = fmap (fmap reads) indexLines :: [[[(Int,String)]]]
                          indices = fmap (fmap $ fst . head) indexReads
                          [maxRow, maxCol] = if nonzeros == 0 then [1,1] else
                                             fmap maximum [fmap (!! 0) indices, fmap (!! 1) indices]
                          valueLines = concat . fmap snd $ splitEntries
                          valueReads = fmap reads valueLines :: [[(Double,String)]]
                          values = fmap (fst . head) valueReads
                          ijxs = zipWith (\[i,j] x -> (i,j,x)) indices values