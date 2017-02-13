// Computing the prestige vector of a social network by power iteration
// Finding the principal eigenvector (the eigenvector corresponding to the largest eigenvalue)
// See http://www.cs.ccsu.edu/~markov/ccsu_courses/dmw2.ppt

public class Eigenvector
{
    static final double epsilon = 0.01;     // convergency threshold
    static double lambda;                   // eigenvalue
    static private int[][] a = { {0,0,1},   // transposed adjaceny matrix
                                 {1,0,0},   // of the social network graph
                                 {1,1,0}, 
                               };
    static double[] p = {1,1,1};            // initial eigenvector

    public static void main (String[] args)
    {
        double[] q;

        System.out.println();
        for (int i=0;i<p.length;i++) System.out.print("p["+i+"]\t\t");
        System.out.println("lambda");
        
        do
        {
            for (int i=0;i<p.length;i++) System.out.format("%f \t",p[i]);  // show the current vector
            System.out.format("%f \n",lambda);                             // and the eigenvalue (lambda)

            q = p;            
            p = AxP(a,q);
            lambda = norm(p);
            p = PxL(p,1/lambda);

        } while (norm(PminusQ(p,q))>epsilon);
   }

// Computes P = A x P (A - matrix, P - column vector)
    public static double[] AxP(int[][] a, double[] p)
    {
        double[] q = new double[p.length];
        double s;
        for (int i=0;i<p.length;i++)
        {
            q[i] = 0;
            for (int j=0;j<p.length;j++)
                q[i] = q[i] + a[i][j]*p[j];
        }
        return q;
    }
    
// Computes P = P x L (P - vector, L - scalar)    
    public static double[] PxL (double[] p, double lambda)
    {
        double[] q = new double[p.length];
        for (int i=0;i<p.length;i++)
            q[i] = p[i]*lambda;
        return q;
    }
    
// Computes P-Q (P and Q - vectors)    
    public static double[] PminusQ (double[] p, double[] q)
    {
        double[] r = new double[p.length];
        for (int i=0;i<p.length;i++)
            r[i] = p[i] - q[i];
        return r;
    }

// Computes Euclidean norm of P
    public static double norm (double[] p)
    {
        double s = 0;
        for (int i=0;i<p.length;i++)
            s = s + p[i]*p[i];
        return Math.sqrt(s);
    }
}
