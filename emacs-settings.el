;; Copy and paste the commands into your ~/.emacs file.

;; Some usefile modes.
(add-hook 'f90-mode-hook 'auto-complete-mode)
(add-hook 'f90-mode-hook 'auto-revert-mode)

;; Add proper comments for f90 mode.
(add-hook 'f90-mode-hook
	  (lambda ()
	    (setq comment-start "!> ")
	    (setq comment-continue "!! ")))

;; Some style.
(setq f90-auto-keyword-case 'downcase-word)

;; We don't like trailing whitespace.
(setq-default before-save-hook 'delete-trailing-whitespace)

;; Electric pair mode...
(electric-pair-mode)

;; Show matching parentheses.
(show-paren-mode)

;; Show trailing whitespace.
(setq-default show-trailing-whitespace t)

;; Show function name in status line.
(which-function-mode)
