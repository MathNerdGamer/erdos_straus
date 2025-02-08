#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <NTL/ZZ.h>
#include <omp.h>

using namespace std::chrono;

// Utility
std::vector<std::uint32_t> primes;
void build_primes( std::string const &, std::uint64_t );
auto timer( steady_clock::time_point const &,
            steady_clock::time_point const & ) -> std::string;

// Erdos-Straus
auto special_cases( NTL::ZZ const & );
auto continued_fraction( NTL::ZZ, NTL::ZZ );
auto generate_continuants( std::vector<std::int64_t> const & );
void erdos_straus( NTL::ZZ const & );
void erdos_straus_fallback( NTL::ZZ const & );

std::array<int, 14> special{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
enum
{
    ZERO_MOD2,
    THREE_MOD4,
    FIVE_MOD8,
    ZERO_MOD3,
    TWO_MOD3,
    ZERO_MOD7,
    THREE_MOD7,
    FIVE_MOD7,
    SIX_MOD7,
    ZERO_MOD5,
    ONE_MOD5,
    TWO_MOD5,
    THREE_MOD5,
    FALLBACK
};

// Config
constexpr bool VERBOSE{ true };
constexpr bool SINGLE{ true };
constexpr bool PRIMES{ false };

auto main( int argc, char **argv ) -> int
{
    NTL::ZZ N{ 1000000 };

    if( argc > 1 )
    {
        std::istringstream input_stream( argv[1] );
        input_stream >> N;

        if( N < 2 )
        {
            std::cout << "N must be at least 2.\n";
            N = NTL::ZZ{ 2 };
        }
    }

    if constexpr( PRIMES )
    {
        std::cout << "Building vector of primes.\n";
        auto const build_start{ high_resolution_clock::now() };
        build_primes( "primes.bin", NTL::conv<unsigned>( N ) );
        auto const build_end{ high_resolution_clock::now() };
        std::cout << "Prime vector of size " << primes.size() << " took "
                  << timer( build_start, build_end ) << " to build.\n";

        std::cout << "Running erdos_straus on each prime up to "
                  << primes.back() << ".\n\n";
    }
    else if constexpr( !SINGLE )
    {
        std::cout << "Running erdos_straus on all integers from 2 to " << N << ".\n\n";
    }

    auto const calculation_start{ high_resolution_clock::now() };

    if constexpr( PRIMES )
    {
        int million = 0;
        std::for_each( std::begin( primes ), std::end( primes ),
                       [](std::int32_t n) { erdos_straus( NTL::ZZ{ n } ); } );
    }
    else if constexpr( SINGLE )
    {
        erdos_straus( N );
    }
    else
    {
        for( NTL::ZZ n{ 2 }; n <= N; ++n )
        {
            erdos_straus( n );
        }
    }

    auto const calculation_end{ high_resolution_clock::now() };

    std::cout << "Finished in " << timer( calculation_start, calculation_end ) << "\n\n";

    if constexpr( !SINGLE )
    {
        std::cout << "In each case below, the numbers are those which were not caught "
                     "in the cases above it.\n";
        std::cout << "Total number = 0 (mod 2): " << special[ZERO_MOD2]  << ".\n";
        std::cout << "Total number = 3 (mod 4): " << special[THREE_MOD4] << ".\n";
        std::cout << "Total number = 5 (mod 8): " << special[FIVE_MOD8]  << ".\n";
        std::cout << "Total number = 0 (mod 3): " << special[ZERO_MOD3]  << ".\n";
        std::cout << "Total number = 2 (mod 3): " << special[TWO_MOD3]   << ".\n";
        std::cout << "Total number = 0 (mod 7): " << special[ZERO_MOD7]  << ".\n";
        std::cout << "Total number = 3 (mod 7): " << special[THREE_MOD7] << ".\n";
        std::cout << "Total number = 5 (mod 7): " << special[FIVE_MOD7]  << ".\n";
        std::cout << "Total number = 6 (mod 7): " << special[SIX_MOD7]   << ".\n";
        std::cout << "Total number = 0 (mod 5): " << special[ZERO_MOD5]  << ".\n";
        std::cout << "Total number = 1 (mod 5): " << special[ONE_MOD5]   << ".\n";
        std::cout << "Total number = 2 (mod 5): " << special[TWO_MOD5] << ".\n";
        std::cout << "Total number = 3 (mod 5): " << special[THREE_MOD5] << ".\n";

        auto total{ std::accumulate( std::begin( special ), std::end( special ), 0 ) - special[FALLBACK] };
        std::string type{ "primes" };

        if constexpr( !PRIMES )
        {
            type = "integers";
        }

        std::cout << "Total number of special cases: " << total << ".\n";
        std::cout << "Total number of " << type << " which had to run through the main algorithm: "
                                                << N - 1 - total << ".\n";
        std::cout << "Total number which used the fallback: " << special[FALLBACK] << ".\n";
    }

    return EXIT_SUCCESS;
}

// Utility Functions
void build_primes( std::string const &file_name, std::uint64_t num_primes )
{   // Builds a vector of primes from a text file.
    std::ifstream in_file{ file_name, std::ios::binary };

    in_file.seekg( 0, std::ios::end );
    std::uint64_t const prime_count{ in_file.tellg() / sizeof( std::uint32_t ) };

    num_primes = std::min( prime_count, num_primes );

    primes.resize( num_primes );
    in_file.seekg(0, std::ios::beg);

    in_file.read( reinterpret_cast<char *>( primes.data() ), num_primes * sizeof( std::uint32_t ) );
}

auto timer( steady_clock::time_point const &start,
            steady_clock::time_point const &end ) -> std::string
{
    auto duration = end - start;

    auto h = duration_cast<hours>( duration );
    duration -= h;

    auto m = duration_cast<minutes>( duration );
    duration -= m;

    auto s = duration_cast<seconds>( duration );
    duration -= s;

    auto ms = duration_cast<milliseconds>( duration );

    std::stringstream elapsed_time;
    bool display = false;

    if( h.count() > 0 )
    {
        elapsed_time << h.count() << "h";
        display = true;
    }

    if( m.count() > 0 || display )
    {
        if( h.count() > 0 && m.count() < 10 )
        {
            elapsed_time << "0";
        }
        elapsed_time << m.count() << "m";
        display = true;
    }

    if( s.count() > 0 || display )
    {
        if( m.count() > 0 && s.count() < 10 )
        {
            elapsed_time << "0";
        }
        elapsed_time << s.count() << "s";
    }

    if( s.count() > 0 && ms.count() < 100 )
    {
        elapsed_time << "0";
    }

    if( s.count() > 0 && ms.count() < 10 )
    {
        elapsed_time << "0";
    }

    elapsed_time << ms.count() << "ms";

    return elapsed_time.str();
}


// Erdos–Straus Functions
auto special_cases( NTL::ZZ const &N )
{   // Handling special cases which were derived in the blog post.
    // Note that this order is necessary, as deriving each formula
    // was conditioned on the previous special case checks failing.
    if( N % 2 == 0 )
    {
        ++special[ZERO_MOD2];

        if constexpr( VERBOSE )
        {
            std::cout << "Special Case: Even.\n";
            std::cout << "4/" << N << " = 1/" << N / 2
                                   << " + 1/" << N
                                   << " + 1/" << N << ".\n";
        }

        return true;
    }

    if( N % 4 == 3 )
    {
        ++special[THREE_MOD4];

        if constexpr( VERBOSE )
        {
            std::cout << "Special Case: 3 (mod 4).\n";
            std::cout << "4/" << N << " = 1/" << ( N + 1 ) / 4
                                   << " + 1/" << N * ( N + 1 ) / 2
                                   << " + 1/" << N * ( N + 1 ) / 2 << ".\n";
        }

        return true;
    }

    if( N % 8 == 5 )
    {
        ++special[FIVE_MOD8];
        if constexpr( VERBOSE )
        {
            std::cout << "Special Case: 5 (mod 8).\n";
            std::cout << "4/" << N << " = 1/" << ( N + 3 ) / 4
                                   << " + 1/" << N * ( ( N + 3 ) / 8 )
                                   << " + 1/" << N * ( ( N + 3 ) / 4 ) << ".\n";
        }

        return true;
    }

    switch( N % 3 )
    {
        case 0:
        {
            ++special[ZERO_MOD3];

            if constexpr( VERBOSE )
            {
                std::cout << "Special Case: Divisible by 3.\n";
                std::cout << "4/" << N << " = 1/" << N / 3
                                       << " + 1/" << 2 * N
                                       << " + 1/" << 2 * N << ".\n";
            }

            return true;
        }

        case 2:
        {
            ++special[TWO_MOD3];

            if constexpr( VERBOSE )
            {
                std::cout << "Special Case: 2 (mod 3).\n";
                std::cout << "4/" << N << " = 1/" << ( N + 1 ) / 3
                                       << " + 1/" << N
                                       << " + 1/" << N * ( N + 1 ) / 3 << ".\n";
            }

            return true;
        }
    }

    switch( N % 7 )
    {
        case 0:
        {
            ++special[ZERO_MOD7];

            if constexpr( VERBOSE )
            {
                std::cout << "Special Case: Divisible by 7.\n";
                std::cout << "4/" << N << " = 1/" << 2 * ( N / 7 )
                                       << " + 1/" << 28 * ( N / 7 )
                                       << " + 1/" << 28 * ( N / 7 ) << ".\n";
            }

            return true;
        }

        case 3:
        {
            ++special[THREE_MOD7];

            if constexpr( VERBOSE )
            {
                std::cout << "Special Case: 3 (mod 7).\n";
                std::cout << "4/" << N << " = 1/" << ( 2 * N + 1 ) / 7
                                       << " + 1/" << 2 * N
                                       << " + 1/" << N * ( 2 * N + 1 ) / 7 << ".\n";
            }

            return true;
        }

        case 5:
        {
            ++special[FIVE_MOD7];

            if constexpr( VERBOSE )
            {
                std::cout << "Special Case: 5 (mod 7).\n";
                std::cout << "4/" << N << " = 1/" << 2 * ( N + 2 ) / 7
                                       << " + 1/" << 2 * N
                                       << " + 1/" << N * ( N + 2 ) / 7 << ".\n";
            }

            return true;
        }

        case 6:
        {
            ++special[SIX_MOD7];

            if constexpr( VERBOSE )
            {
                std::cout << "Special Case: 6 (mod 7).\n";
                std::cout << "4/" << N << " = 1/" << 2 * ( N + 1 ) / 7
                                       << " + 1/" << 2 * N
                                       << " + 1/" << 2 * N * ( N + 1 ) / 7 << ".\n";
            }

            return true;
        }
    }

    switch( N % 5 )
    {
        case 0:
        {
            ++special[ZERO_MOD5];

            if constexpr( VERBOSE )
            {
                std::cout << "Special Case: Divisible by 5.\n";
                std::cout << "4/" << N << " = 1/" << 2 * ( N / 5 )
                                       << " + 1/" << N
                                       << " + 1/" << 2 * N << ".\n";
            }

            return true;
        }

        case 1:
        {
            ++special[ONE_MOD5];

            if constexpr( VERBOSE )
            {
                std::cout << "Special Case: 1 (mod 5).\n";
                std::cout << "4/" << N << " = 1/" << 4 * ( N + 4 ) / 15
                                       << " + 1/" << 4 * N
                                       << " + 1/" << N * ( N + 4 ) / 15 << ".\n";
            }

            return true;
        }

        case 2:
        {
            ++special[TWO_MOD5];

            if constexpr( VERBOSE )
            {
                std::cout << "Special Case: 2 (mod 5).\n";
                std::cout << "4/" << N << " = 1/" << 2 * ( 2 * N + 1 ) / 15
                                       << " + 1/" << 4 * N
                                       << " + 1/" << 4 * N * ( 2 * N + 1 ) / 15 << ".\n";
            }

            return true;
        }

        case 3:
        {
            ++special[THREE_MOD5];

            if constexpr( VERBOSE )
            {
                std::cout << "Special Case: 3 (mod 5).\n";
                std::cout << "4/" << N << " = 1/" << 4 * ( N + 2 ) / 15
                                       << " + 1/" << 4 * N
                                       << " + 1/" << 2 * N * ( N + 2 ) / 15 << ".\n";
            }

            return true;
        }
    }

    // Anything past this point must be handled by the full algorithm.
    if constexpr( VERBOSE )
    {
        std::cout << "Not among the special cases covered; full algorithm will be performed.\n";
    }
    return false;
}

auto continued_fraction( NTL::ZZ numerator, NTL::ZZ denominator )
{   // Continued fraction a[0] + 1/(a[1] + 1/(a[2] + ...))
    // Stored more succintly as [ a[0], a[1], a[2], ... ]
    std::vector<std::int64_t> coefficients;

    while( denominator != 0 )
    {
        auto quotient{ numerator / denominator }; // a[k]
        coefficients.push_back( NTL::conv<int>( quotient ) );

        auto remainder{ numerator % denominator }; // Save remainder of numerator / denominator in temp variable
        numerator = denominator;
        denominator = remainder;
    }
    return coefficients;
}

auto generate_continuants( std::vector<std::int64_t> const &continued_fraction )
{   // We only need the denominators for our purposes.
    std::vector<NTL::ZZ> continuants;

    NTL::ZZ A_prev{ 0 }, A_curr{ 1 }, // A[-1] = 0, A[0] = 1
            B_prev{ 1 }, B_curr{ 0 }; // B[-1] = 1, B[0] = 0

    for( int coefficient : continued_fraction )
    {   // A[k + 1] = a[k + 1] * A[k] + A[k - 1]
        NTL::ZZ A_next = A_prev + NTL::ZZ{ coefficient } *A_curr;

        // B[k + 1] = a[k + 1] * B[k] + B[k - 1]
        NTL::ZZ B_next = B_prev + NTL::ZZ{ coefficient } *B_curr;

        // x[k] = A[k] / B[k]
        // Save the denominator
        continuants.emplace_back( B_next );

        // k -> k+1
        A_prev = A_curr;
        A_curr = A_next;
        B_prev = B_curr;
        B_curr = B_next;
    }

    return continuants;
}

void erdos_straus( NTL::ZZ const &N )
{   // Finding x,y,z such that 4/N = 1/x + 1/y + 1/z
    if( N < 2 )
    {
        return;
    }

    if( special_cases(N) )
    {
        return;
    }

    auto const x_min{ ( N + 3 ) / 4 };
    // ceil(N/4) = floor((N+3) / 4)
    auto const x_max{ ( 3 * N ) / 4 };
    // floor(3N/4)

    bool found{ false };
    auto thread_count = omp_get_max_threads();

    #pragma omp parallel shared(found)
    {
        auto const thread_num = omp_get_thread_num();
        for( auto x{ x_min + thread_num }; x <= x_max && !found; x += thread_count )
        {   // Considering 4/N = 1/x + 1/y + 1/z,
            // we have 4/N - 1/x     =
            //         (4x - N)/(Nx) = 1/y + 1/z
            auto const numerator{ 4 * x - N };
            auto const denominator{ N * x };
            if( numerator <= 0 )
            {
                continue;
            }

            // Generate the denominator continuants for the convergents of the continued fraction
            auto const cf{ continued_fraction( numerator, denominator ) };
            auto const continuants{ generate_continuants( cf ) };

            for( auto const &B : continuants )
            {   // y ~= denominator of a convergent
                // (1/y should be close to (4x - N)/(Nx))
                if( found )
                {
                    break;
                }

                auto const y{ B };
                if( y < x || y >( 2 * N * x ) / numerator )
                {   // If y is too small or big, skip it.
                    continue;
                }

                // 1/z = 4/N - 1/x - 1/y = numerator/denominator - 1/y
                //     = (y * numerator - denominator) / (y * denominator)
                auto const z_numerator{ denominator * y };
                auto const z_denominator{ numerator * y - denominator };

                if( z_denominator <= 0 || z_numerator % z_denominator != 0 )
                {
                    continue;
                }

                auto const z{ z_numerator / z_denominator };

                if( z >= y )
                {
                    #pragma omp critical
                    if( !found )
                    {
                        found = true;
                        if constexpr( VERBOSE )
                        {
                            std::cout << "4/" << N << " = 1/" << x
                                                   << " + 1/" << y
                                                   << " + 1/" << z << ".\n";
                        }
                    }
                    break;
                }
            }
        }
    }

    if( !found )
    {
        if constexpr( VERBOSE )
        {
            std::cout << "Entering the fallback.\n";
        }
        erdos_straus_fallback( N );
    }

    return;
}

void erdos_straus_fallback( NTL::ZZ const &N )
{   // Finding x,y,z such that 4/N = 1/x + 1/y + 1/z
    ++special[FALLBACK];

    std::cout << "Entered fallback at N = " << N << ".\n";
    return;

    auto const x_min{ ( N + 3 ) / 4 };
    // ceil(N/4) = floor((N+3) / 4)
    auto const x_max{ ( 3 * N ) / 4 };
    // floor(3N/4)
    
    bool found{ false };
    auto thread_count = omp_get_max_threads();
    
    #pragma omp parallel shared(found)
    {
        auto const thread_num = omp_get_thread_num();
        for( auto x{ x_min + thread_num }; x <= x_max && !found; x += thread_count )
        {   // ceil(N/4) <= x <= floor(3N/4)
            auto const y_denominator{ 4 * x - N }; // 4x - N
            if( y_denominator <= 0 )
            {
                continue;
            }
    
            auto const y_max{ ( 2*N*x ) / y_denominator };
            // floor(2Nx / (4x-N))
    
            for( auto y{ x }; y <= y_max && !found; ++y )
            {   // x <= y <= floor(2Nx / (4x-N))
                auto const z_denominator{ 4 * x * y - N * ( x + y ) };
                // 4xy - N(x+y)
                if( z_denominator <= 0 )
                {
                    continue;
                }
    
                auto const z_numerator{ N * x * y };
                // Nxy
    
                if( z_numerator % z_denominator != 0 )
                {   // z must be an integer.
                    continue;
                }
    
                auto const z{ z_numerator / z_denominator };
                // z = Nxy / (4xy - N(x+y))
    
                if( z >= y )
                {
                    #pragma omp critical
                    if( !found )
                    {
                        found = true;
                        if constexpr( VERBOSE )
                        {
                            std::cout << "4/" << N << " = 1/" << x
                                                   << " + 1/" << y
                                                   << " + 1/" << z << ".\n";
                        }
                    }
                    break;
                }
            }
        }
    }
    
    
    if( !found )
    {
        std::cout << "Counterexample at N = " << N << "?\n";
    }
    
    return;
}
