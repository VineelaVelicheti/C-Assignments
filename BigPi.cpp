#include <iostream>
#include <iomanip>
#include <mpir.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>

using namespace std;

const int MAX_ITERATIONS = 100;
const int PLACES         = 1000;        // desired decimal places
const int PRECISION      = PLACES + 1;  // +1 for the digit 3 before the decimal

const int BASE       = 10;  // base 10 numbers
const int BIT_COUNT  = 8;   // bits per machine word

const int BLOCK_SIZE = 10;                // print digits in blocks
const int LINE_SIZE  = 100;               // digits to print per line
const int LINE_COUNT = PLACES/LINE_SIZE;  // lines to print
const int GROUP_SIZE = 5;                 // line grouping size

/**
 * Compute the cube root of a positive integer.
 * @param x where to store the result.
 * @param a the number whose cube root to compute.
 */
void cube_root(mpf_t& x, const mpf_t a);
/* Compute the value of pi using nonic algorithm
 * @param  stores the value of pi after computation*/
void pi_nonic_decimals(mpf_t& x);
/* Prints the value of pi in required format to a file*/
void print_pi_decimals(char * pi_value,mpf_t pi);
/***** Add more functions as necessary. *****/

/**
 * The main.
 */
int main()
{
    mpf_set_default_prec(BIT_COUNT*PRECISION);  // precision in bits
    mpf_t pi; mpf_init(pi);
    pi_nonic_decimals(pi);
    char * pi_value;
    mp_exp_t exp;
    pi_value=mpf_get_str(NULL,&exp,BASE,PRECISION,pi);
    print_pi_decimals(pi_value,pi);
    return 0;
}

void pi_nonic_decimals(mpf_t& x)
{
	mpf_t a,r,s,t,u,v,w; //Declare all Nonic variables needed for implementing Nonic algorithm
	mpf_t a01,a02,a03,v01,s01; //Temporary variables needed for calculating the actual Nonic Values
	// Initialize all the variables
	mpf_init(a);
	mpf_init(r);
	mpf_init(s);
	mpf_init(t);
	mpf_init(u);
	mpf_init(v);
	mpf_init(w);
	mpf_init(v01);
	mpf_init(a01);
	mpf_init(a02);
	mpf_init(a03);
	mpf_init(s01);
	//Nonic Algorithm Implementation
	// Calculating the initial values for a,r,s
	mpf_set_ui(a,1);
	mpf_div_ui(a,a,3);
	mpf_set_ui(r,3);
	mpf_sqrt (r,r);
	mpf_sub_ui(r,r,1);
	mpf_div_ui(r,r,2);
	mpf_pow_ui(s,r,3);
	mpf_ui_sub(s,1,s);
    cube_root(s,s);
    int n=0;
    //loop for a finite iterations based on the decimal points needed
    while(n<MAX_ITERATIONS)
    {
    		//Calculate value for t
    		mpf_mul_ui(t,r,2);
    		mpf_add_ui(t,t,1);
    		//Calculate value for u
    		mpf_pow_ui(u,r,2);
    		mpf_add(u,u,r);
    		mpf_add_ui(u,u,1);
    		mpf_mul_ui(u,u,9);
    		mpf_mul(u,u,r);
    		cube_root(u,u);
    		//Calculate value for v
    		mpf_pow_ui(v,u,2);
    		mpf_mul(v01,t,u);
    		mpf_add(v,v,v01);
    		mpf_pow_ui(v01,t,2);
        mpf_add(v,v,v01);
        //Calculate value for w
        mpf_pow_ui(w,s,2);
        mpf_add(w,w,s);
        mpf_add_ui(w,w,1);
        mpf_mul_ui(w,w,27);
        mpf_div(w,w,v);
        //Update value of a
        mpf_mul(a01,w,a);
        mpf_ui_sub(a02,1,w);
        mpf_set_ui(a03,3);
        int k=2*n;
        if(k==0)
        {
        	mpf_ui_div(a03,1,a03);
        }
        else
        {
        	mpf_pow_ui(a03,a03,k-1);
        }
        mpf_mul(a02,a02,a03);
        mpf_add(a01,a01,a02);
        mpf_set(a,a01);
        // Update initial value of s
        mpf_ui_sub(s,1,r);
        mpf_pow_ui(s,s,3);
        mpf_mul_ui(s01,u,2);
        mpf_add(s01,s01,t);
        mpf_mul(s01,s01,v);
        mpf_div(s,s,s01);
        //Update initial value of r
        mpf_pow_ui(r,s,3);
        mpf_ui_sub(r,1,r);
        cube_root(r,r);
        n++;
    }
    // Update the value of pointer x so that the value in pi in updated
    mpf_ui_div(x,1,a);
    //clear all values
    	mpf_clear(a);
    	mpf_clear(r);
    	mpf_clear(s);
    	mpf_clear(t);
    	mpf_clear(u);
    	mpf_clear(v);
    	mpf_clear(w);
    	mpf_clear(a01);
    	mpf_clear(a02);
    	mpf_clear(a03);
    	mpf_clear(v01);
    	mpf_clear(s01);
}

void cube_root(mpf_t& x, const mpf_t a)
{
    mpf_t cube_root_initial_approx;
	mpf_t x_value,x_prev_value;
	mpf_t x1,x2,x3,a1;
    // Initialize all the values to 0
    mpf_init(cube_root_initial_approx);
    mpf_init(x_value);
    mpf_init(x_prev_value);
	mpf_init(x1);
	mpf_init(x2);
	mpf_init(x3);
	mpf_init(a1);
	mpf_div_ui(cube_root_initial_approx,a,3);// Setting the initial estimate value to be ONE_THIRD of the value for better performance
	mpf_set(x_value, cube_root_initial_approx);
	int i=0;
	while(mpf_cmp(x_value,x_prev_value)!=0 && i<MAX_ITERATIONS) // loop till the max number of iterations or till the previous value and the current value are same
	{
		mpf_pow_ui(x1,x_value,3);
		mpf_mul_ui(a1,a,2);
		mpf_add(x2,x1,a1);
		mpf_mul_ui(x1,x1,2);
		mpf_add(x3,x1,a);
		mpf_div(x2,x2,x3);
		mpf_set(x_prev_value,x_value);
		mpf_mul(x_value,x_value,x2);
		i++;
	}
	mpf_set(x,x_value);
	//clear all values
	mpf_clear(x_value);
	mpf_clear(x_prev_value);
	mpf_clear(cube_root_initial_approx);
	mpf_clear(x1);
	mpf_clear(x2);
	mpf_clear(x3);
	mpf_clear(a1);
}

void print_pi_decimals(char * pi_value,mpf_t pi)
{
	mpir_si d=mpf_get_si(pi);
	ofstream outstream;
    outstream.open("BigpiOutput.txt");
    outstream<<d<<".";
    pi_value=pi_value+1;
    int ls=0,gs=0,lc=0;
	while(lc<LINE_COUNT)
	{
	    	ls=0;
		while(ls<LINE_SIZE)
		{
			for(int i=0;i<BLOCK_SIZE;i++)
			{
			  outstream<<pi_value[i];
			  ls++;
			}
			pi_value=pi_value+BLOCK_SIZE;
			outstream<<" ";
		}
		outstream<<endl;
		outstream<<"  ";
		gs++;
		if(gs==GROUP_SIZE)
		{
			outstream<<endl;
			outstream<<"  ";
		}
		lc++;
	}
	outstream.close();

}
