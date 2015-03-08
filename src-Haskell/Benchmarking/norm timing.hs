import MatrixTree
import SpAMM
import System.IO (hClose, hPutStr, openFile, IOMode(WriteMode))
import System.Process
import System.TimeIt

main = do let sizes = fmap (2^) [11..20]
          putStrLn "norm"
          mapM_ (doTiming setNorm) sizes

doTiming :: (MatrixTree -> Norm) -> Int -> IO ()
doTiming op size = do let filePath = show size ++ ".txt"
                      createMatrix filePath size
                      tree <- readTreeFromMatrixMarket filePath
                      (time, norm) <- timeItT . return $ op tree
                      putStrLn $ show size ++ " " ++ show time
                      callProcess "rm" [filePath]

createMatrix :: FilePath -> Int -> IO ()
createMatrix filePath size = do
             handle <- openFile filePath WriteMode
             matrix <- readProcess "python" ["../generate-matrix.py", "-N " ++ show size] []
             hPutStr handle matrix
             hClose handle
