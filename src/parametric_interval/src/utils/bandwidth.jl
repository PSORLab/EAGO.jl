# Convert Sparse Matrix to Band Matrix
function sparse_bandwidth(A::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti<:Integer}
     # get start column
     kl::Int = 1
     ku::Int = 1
     colind::Int = 1
     for i=1:(length(A.colptr)-1)
          for j=(A.colptr[i]):(A.colptr[i+1]-1)
               del = colind-A.rowval[j]
               if del>0
                    (del>ku) && (ku = del)
               else
                    (del<-kl) && (kl = -del)
               end
          end
          colind +=1
     end
     return kl,ku
end
