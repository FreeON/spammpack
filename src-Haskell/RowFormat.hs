module RowFormat
( MatrixTree
, readTreeFromRowFormat
, writeTreeToRowFormat
) where

-- reads/writes a MatrixTree to a file format that lists rows in order on separate lines
-- with entries on each row separated by spaces

import Data.List (intersperse)
import MatrixTree (MatrixTree, RowList, rowListToTree, treeToRowList)
import System.IO (hClose, hPutStrLn, openFile, IOMode(WriteMode))

readTreeFromRowFormat :: FilePath -> IO MatrixTree
readTreeFromRowFormat filePath = readRowListFromRowFormat filePath >>= (return . rowListToTree)

writeTreeToRowFormat :: MatrixTree -> FilePath -> IO ()
writeTreeToRowFormat tree filePath = writeRowListToRowFormat (treeToRowList tree) filePath

readRowListFromRowFormat :: FilePath -> IO RowList
readRowListFromRowFormat filePath = do contents <- readFile filePath
                                       let rowList = fmap (fmap read . words) (lines contents)
                                       if isValidRowFormat rowList then return rowList
                                       else error $ "Row format file " ++ filePath ++ " is invalid"

writeRowListToRowFormat :: RowList -> FilePath -> IO ()
writeRowListToRowFormat rowList filePath =
                        do handle <- openFile filePath WriteMode
                           mapM_ (hPutStrLn handle) $ fmap (combineStrings . (fmap show)) rowList
                           hClose handle

isValidRowFormat :: RowList -> Bool
isValidRowFormat rows = not (null rows) && not (null $ head rows) && all sameLength rows
                        where sameLength row = length row == length (head rows)

combineStrings :: [String] -> String
combineStrings = concat . (intersperse " ")