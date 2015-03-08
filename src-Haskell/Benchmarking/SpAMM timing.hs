import MatrixTree
import SpAMM
import System.IO (hClose, hPutStrLn, openFile, IOMode(WriteMode))
import System.Process
import System.TimeIt

main = do let sizes = fmap (2^) [11..20]
          putStrLn "treeAdd"
          mapM_ (doTiming treeAdd) sizes
          putStrLn "treeMult"
          mapM_ (doTiming treeMult) sizes

doTiming :: (MatrixTree -> MatrixTree -> MatrixTree) -> Int -> IO ()
doTiming op size = do let filePath1 = "fst" ++ show size ++ ".txt"
                      let filePath2 = "snd" ++ show size ++ ".txt"
                      createMatrix filePath1 size
                      createMatrix filePath2 size
                      fstTree <- readTreeFromMatrixMarket filePath1
                      sndTree <- readTreeFromMatrixMarket filePath2
                      (time, tree) <- timeItT . return $ op fstTree sndTree
                      putStrLn $ show size ++ " " ++ show time
                      callProcess "rm" [filePath1, filePath2]

createMatrix :: FilePath -> Int -> IO ()
createMatrix filePath size = do
             handle <- openFile filePath WriteMode
             matrix <- readProcess "python" ["../generate-matrix.py", "-N " ++ show size] []
             hPutStrLn handle matrix
             hClose handle
