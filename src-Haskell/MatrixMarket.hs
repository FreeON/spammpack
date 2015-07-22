module MatrixMarket
( MatrixList
, mmReadFile
, mmWriteFile
) where

-- reads/writes MatrixMarket files to a generic matrix data type

import Data.Char (toLower)
import Data.Complex
import Data.List (intersperse, nub)
import qualified Data.Map as Map (fromList, lookup, Map)
import Data.Maybe (fromJust, isNothing)
import MatrixList (MatrixList)
import System.IO (hClose, hPutStrLn, openFile, IOMode(WriteMode))

data MMFormat = Array | Coordinate
data MMType = MMInt | MMReal | MMComplex

class MMReadable a where
      mmRead :: String -> a
      mmReads :: String -> [(a, String)]

instance MMReadable Int where
         mmRead = read
         mmReads = reads

instance MMReadable Float where
         mmRead = read
         mmReads = reads

instance MMReadable Double where
         mmRead = read
         mmReads = reads

instance (MMReadable a, RealFloat a) => MMReadable (Complex a) where
         mmRead strz = a :+ b where [a, b] = fmap mmRead (words strz)

         mmReads str | len == 2  = if null xread || null yread then []
                                   else [(x :+ y, "")]
                     | otherwise = []
                     where len = length (words str)
                           [xread, yread] = fmap mmReads (words str)
                           [x, y] = fmap (fst . head) [xread, yread]

instance MMReadable MMFormat where
         mmReads name =
                 case map toLower name of "array"      -> [(Array, "")]
                                          "coordinate" -> [(Coordinate, "")]
                                          _            -> []

         mmRead = fst . head . mmReads

instance MMReadable MMType where
         mmReads name =
                 case map toLower name of "integer" -> [(MMInt, "")]
                                          "real"    -> [(MMReal, "")]
                                          "complex" -> [(MMComplex, "")]
                                          _         -> []

         mmRead = fst . head . mmReads

class MMShowable a where
      mmShow :: a -> String

instance MMShowable Int where
         mmShow = show

instance MMShowable Float where
         mmShow = show

instance MMShowable Double where
         mmShow = show

instance (MMShowable a, RealFloat a) => MMShowable (Complex a) where
         mmShow z = concat . intersperse " " . fmap mmShow $ [realPart z, imagPart z]

instance MMShowable MMFormat where
         mmShow Array = "array"
         mmShow Coordinate = "coordinate"

instance MMShowable MMType where
         mmShow MMInt = "integer"
         mmShow MMReal = "real"
         mmShow MMComplex = "complex"

mmReadFile :: FilePath -> IO (MatrixList Double)
mmReadFile filePath = readFile filePath >>=
                      (return . mmToMatrixList filePath)

mmWriteFile :: MatrixList Double -> String -> FilePath -> IO ()
mmWriteFile matrixList format filePath = do
            let action = Map.lookup format mmWriteCmds
            let list = maybe (error $ "Unrecognized format " ++ format)
                             ($ matrixList) action
            handle <- openFile filePath WriteMode
            hPutStrLn handle $ "%%MatrixMarket matrix " ++ format ++ " real general"
            mapM_ (hPutStrLn handle) list
            hClose handle

mmToMatrixList :: FilePath -> String -> MatrixList Double
mmToMatrixList filePath contents
      | null (words contents) = error $ filePath ++ " is empty"
      | length first < 4      = error $ filePath ++ " header is invalid"
      | null rest             = error $ filePath ++ " has no contents"
      | otherwise             = maybe errormsg (($ filePath) . ($ rest)) action
      where fileLines = lines contents
            first = words . head $ fileLines
            rest = filter ((/='%') . head . head) . filter (not . null) .
                   fmap words $ tail fileLines
            action = Map.lookup (fmap toLower (first !! 2)) mmReadCmds
            errormsg = error $ filePath ++ " has unrecognized format " ++
                               (first !! 2)

mmReadCmds :: Map.Map String ([[String]] -> FilePath -> MatrixList Double)
mmReadCmds = Map.fromList [
               ("array", arrayToMatrixList)
             , ("coordinate", coordinateToMatrixList)
             ]

arrayToMatrixList :: [[String]] -> FilePath -> MatrixList Double
arrayToMatrixList lineList filePath
     | length firstLine /= 2             = error $ filePath ++ " matrix size invalid"
     | any ((/= 1) . length) valueLines  = error $ filePath ++ " has invalid value lines"
     | any null size                     = error $ filePath ++ " matrix size invalid"
     | length valueList /= nrows * ncols = error $ filePath ++ " has wrong matrix size"
     | any null valueReads               = error $ filePath ++ " has invalid values"
     | otherwise                         = (nrows, ncols, ijxs)
     where (firstLine, valueLines) = (head lineList, tail lineList)
           size = fmap reads firstLine :: [[(Int, String)]]
           [nrows, ncols] = fmap (fst . head) size
           valueList = concat valueLines
           valueReads = fmap reads valueList :: [[(Double, String)]]
           values = fmap (fst . head) valueReads
           indices = [(i, j) | j <- [1..ncols], i <- [1..nrows]]
           ijxs = fmap (\((i, j), x) -> (i, j, x)) . filter ((/= 0) . snd) .
                  zip indices $ values

coordinateToMatrixList :: [[String]] -> FilePath -> MatrixList Double
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
           firstLineNums = fmap reads firstLine :: [[(Int, String)]]
           firstLineVals = fmap (fst . head) firstLineNums
           [nrows, ncols] = take 2 firstLineVals
           nonzeros = firstLineVals !! 2
           splitEntries = fmap (splitAt 2) entryLines
           indexLines = fmap fst splitEntries
           indexReads = fmap (fmap reads) indexLines :: [[[(Int, String)]]]
           indices = fmap (fmap $ fst . head) indexReads
           [maxRow, maxCol] = if nonzeros == 0 then [1,1] else
                              fmap maximum [fmap (!! 0) indices, fmap (!! 1) indices]
           valueLines = concat . fmap snd $ splitEntries
           valueReads = fmap reads valueLines :: [[(Double, String)]]
           values = fmap (fst . head) valueReads
           ijxs = zipWith (\[i,j] x -> (i,j,x)) indices values

mmWriteCmds :: Map.Map String (MatrixList Double -> [String])
mmWriteCmds = Map.fromList [
                ("array", matrixListToArray)
              , ("coordinate", matrixListToCoordinate)
              ]

matrixListToArray :: MatrixList Double -> [String]
matrixListToArray (m, n, ijxs) = (joinStr [show m, show n]):showVals
                   where pairList = fmap (\(i, j, x) -> ((i,j), x)) ijxs
                         hashTable = Map.fromList pairList
                         indices = [(i, j) | j <- [1..n], i <- [1..m]]
                         showVals = fmap (show . value) indices
                         value pair | isNothing check = 0
                                    | otherwise       = fromJust check
                                    where check = Map.lookup pair hashTable

matrixListToCoordinate :: MatrixList Double -> [String]
matrixListToCoordinate (m, n, ijxs) =
          (joinStr [show m, show n, show $ length ijxs]):(fmap pairToStr ijxs)

joinStr :: [String] -> String
joinStr = concat . intersperse " "

pairToStr :: (Show a, Show b, Show c) => (a,b,c) -> String
pairToStr = joinStr . (\(i,j,x) -> [show i, show j, show x])
