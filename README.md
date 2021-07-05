# mt19937-reversible

This is an educational-, research-, and CTF-oriented implementation of the
Mersenne Twister PRNG, providing features such as recovering internal state from
outputs, rewinding the internal state, or recovering the seed from a freshly
seeded internal state.

There are may implementations of Mersenne Twister on the Internet, and even
may implementations that allow cloning internal state. This one is better
because:

* Both 32-bit MT19937 and 64-bit MT19937-64 are implemented
* A rewind feature is provided to "turn back time" on the PRNG
* The value of the seed can be recovered from a freshly-seeded state
* The code is written in a conceptually simple manner, making it easier to 
reason about. Performance does take a hit, though.

An example that demonstrates cloning a standard library `random.Random` is
provided in the code, and also below:

    def test_with_stdlib_random():
        import random
        stdlib_random = random.Random()
        mt = MT19937()
        output_for_cloning = [stdlib_random.getrandbits(mt.w) for _ in range(mt.n)]
        mt.clone_state_from_output(output_for_cloning)
        for _ in range(10000):
            assert stdlib_random.getrandbits(mt.w) == mt.get_next_random()
        print('Test successful')

Have fun!
