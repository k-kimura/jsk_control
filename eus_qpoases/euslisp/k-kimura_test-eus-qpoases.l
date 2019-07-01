(load "package://eus_qpoases/euslisp/test-eus-qpoases.l")

;; sample for lp
;;max 2x1 + x2
;;    => [2 1] [x1 x2]'
;;s.t -x1 + 2x2 ≤ 4
;;    2x1 - x2 ≤ 6
;;    -x1 - x2 ≤ -2
;;    =>
;;    |-1 2||x1 x2|'<=4
;;    |2 -1||x1 x2|'<=6
;;    |-1 -1||x1 x2|'<=-2
;;=>[5.33333 4.66667]
(defun k-kimura-demo-eus-qpOASES5
  (&key (solve-lp-mode :lp))
  "Demo for LP with inequality."
  ;; (let* ((X0 (float-vector 2.6666 0 0 0))
  (let* ((X0 (float-vector 2.6666 0))
         (q (scale -1 (float-vector 2 1)))
         (A (make-matrix
             3 2
             (list (list -1 2)
                   (list 2 -1)
                   (list -1 -1))))
         (ubA (float-vector 4 6 -2))
         (ret (solve-qpoases-lp
               :initial-state x0
               :eval-coeff-vector q
               ;; :state-min-vector (float-vector 0 0)
               :inequality-matrix A
               :inequality-max-vector ubA
               :solve-lp-mode solve-lp-mode
               )))
    (warn ";; solution ~A => ~A~%"
          ret
          ;; (setq ret
          ;;       (eps-v= (float-vector 22.5 22.5) ret))
          (float-vector 5.33333 4.66667)
          )
    ;; correct solution
    ret
    ))


(defun calc-stability-margin-max-point-from-support-polygon-vertices (spv-list &key (epsilon 0.0)) ;; epsilon>=0
;; (defun make-dist-xy-each-line-list-from-support-polygon-vertices
  (let* (xi yi xi1 yi1 naihou-inequality-matrix naihou-inequality-max-vector di ai bi ci abc-dist-list aj bj cj stability-margin-max-list stability-margin-max-list-final)
    (dotimes (i (length spv-list))
      (setq xi (elt (elt spv-list i) 0))
      (setq yi (elt (elt spv-list i) 1))
      (if (< i (- (length spv-list) 1))
          (progn
            (setq xi1 (elt (elt spv-list (+ i 1)) 0))
            (setq yi1 (elt (elt spv-list (+ i 1)) 1)))
        (progn
          (setq xi1 (elt (car spv-list) 0))
          (setq yi1 (elt (car spv-list) 1)))
        )
      (push (list (- yi1 yi) (- (- xi1 xi))) naihou-inequality-matrix)
      (push (- (* xi yi1) (* xi1 yi)) naihou-inequality-max-vector)
      (setq di (sqrt (+ (expt (- xi1 xi) 2) (expt (- yi1 yi) 2))))
      (setq ai (- (/ (- yi1 yi) di)))
      (setq bi (/ (- xi1 xi) di))
      (setq ci (/ (- (* xi yi1) (* xi1 yi)) di))
      ;;;; dist-xy-line-i = ai*x + bi*y + ci
      (push (float-vector ai bi ci) abc-dist-list)
      )
    (setq abc-dist-list (reverse abc-dist-list))
    (warn "abc-dist-list: ~a~%~%" abc-dist-list)
    (dotimes (i (length abc-dist-list))
      (let* (lp-eval-coeff-vector-i lp-inequality-matrix-i lp-inequality-max-vector-i stability-margin-max-point-i stability-margin-max-point-x-i stability-margin-max-point-y-i stability-margin-max-dist-i)
        (setq ai (elt (elt abc-dist-list i) 0))
        (setq bi (elt (elt abc-dist-list i) 1))
        (setq ci (elt (elt abc-dist-list i) 2))
        (setq lp-eval-coeff-vector-i (float-vector ai bi))
        (warn "lp-eval-coeff-vector-~a = ~a~%" i lp-eval-coeff-vector-i)
        (dotimes (j (length abc-dist-list))
          (unless (= j i)
            (setq aj (elt (elt abc-dist-list j) 0))
            (setq bj (elt (elt abc-dist-list j) 1))
            (setq cj (elt (elt abc-dist-list j) 2))
            (push (list (- ai aj) (- bi bj)) lp-inequality-matrix-i)
            (push (- (- ci cj)) lp-inequality-max-vector-i)
            )
          )
        ;; (setq lp-inequality-matrix-i (append lp-inequality-matrix-i naihou-inequality-matrix))
        ;; (setq lp-inequality-max-vector-i (append lp-inequality-max-vector-i naihou-inequality-max-vector))
        (setq lp-inequality-matrix-i (append lp-inequality-matrix-i naihou-inequality-matrix (list (list 0 -1) (list 0 1))))
        (setq lp-inequality-max-vector-i (append lp-inequality-max-vector-i naihou-inequality-max-vector (list epsilon epsilon)))
        (setq lp-inequality-matrix-i (reverse lp-inequality-matrix-i))
        (setq lp-inequality-max-vector-i (coerce (reverse lp-inequality-max-vector-i) float-vector))
        (warn "lp-inequality-matrix-~a = ~a~%" i lp-inequality-matrix-i)
        (warn "lp-inequality-max-vector-~a = ~a~%" i lp-inequality-max-vector-i)
        (setq stability-margin-max-point-i (solve-qpoases-lp
                                            :initial-state (float-vector 0 0)
                                            :eval-coeff-vector (scale -1 lp-eval-coeff-vector-i)
                                            ;; :state-min-vector (float-vector 0 0)
                                            :inequality-matrix (make-matrix
                                                                (length lp-inequality-matrix-i) 2
                                                                lp-inequality-matrix-i)
                                            :inequality-max-vector lp-inequality-max-vector-i
                                            :solve-lp-mode :lp
                                            ))
        (warn "stability-margin-max-point-~a = ~a~%~%" i stability-margin-max-point-i)
        (unless (equal stability-margin-max-point-i nil)
          (setq stability-margin-max-point-x-i (elt stability-margin-max-point-i 0))
          (setq stability-margin-max-point-y-i (elt stability-margin-max-point-i 1))
          (setq stability-margin-max-dist-i (+ (* ai stability-margin-max-point-x-i) (* bi stability-margin-max-point-y-i) ci)) ;; dist-xy-line-i = ai*x + bi*y + ci
          ;; (send (float-vector stability-margin-max-point-x-i stability-margin-max-point-y-i 70.0) :draw-on :flush t :color #F(1 0 0)) ;; draw
          (push (float-vector stability-margin-max-point-x-i stability-margin-max-point-y-i stability-margin-max-dist-i) stability-margin-max-list)
          )
        ))
    (setq stability-margin-max-list (reverse stability-margin-max-list))
    (warn "stability-margin-max-list: ~a~%~%" stability-margin-max-list)
    (setq stability-margin-max-list-final (output-zmax-fltvec stability-margin-max-list))
    (warn "stability-margin-max-list-final: ~a~%~%" stability-margin-max-list-final)
    ;; (dolist (fltvec stability-margin-max-list-final)
    ;;   (send (float-vector (elt fltvec 0) (elt fltvec 1) 70.0) :draw-on :flush t :color #F(0 1 0)) ;; draw
    ;;   )
    stability-margin-max-list-final
    ))

(defun output-zmax-fltvec (l)
  (let* (zmax outl)
    (setq zmax (elt (apply 'vmax l) 2))
    (dolist (fltvec l)
      (if (= (elt fltvec 2) zmax)
          (push fltvec outl))
      )
    (setq outl (reverse outl))
    outl
    ))


;; ;; 本当の解は#f(29.2893 29.2893)
;; (defun dame-demo-eus-qpOASES ()
;;   (let* ((lp-eval-coeff-vector-0 #f(0.0 1.0))
;;          (lp-inequality-matrix-0 '((0.0 -100.0) (100.0 100.0) (-100.0 0.0) (0.707107 1.70711) (-1.0 1.0)))
;;          (lp-inequality-max-vector-0 #f(0.0 10000.0 0.0 70.7107 0.0))
;;          )
;;     (solve-qpoases-lp
;;      :initial-state (float-vector 0 0)
;;      :eval-coeff-vector (scale -1 lp-eval-coeff-vector-0)
;;      ;; :state-min-vector (float-vector 0 0)
;;      :inequality-matrix (make-matrix
;;                          (length lp-inequality-matrix-0) 2
;;                          lp-inequality-matrix-0)
;;      :inequality-max-vector lp-inequality-max-vector-0
;;      :solve-lp-mode :lp
;;      )
;;     ))
