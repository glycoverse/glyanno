# unreadable caches fall back and failed rebuilds preserve files

    Code
      suppressMessages(build_glyanno_cache(force = TRUE))
    Condition
      Error in `.build_annotation_view()`:
      ! interrupted build

# build and clear respect a concurrent writer

    Code
      clear_glyanno_cache()
    Condition
      Error in `clear_glyanno_cache()`:
      ! Cannot clear the cache while a build lock exists.

# a failed replacement restores the previous cache

    Code
      .replace_cache_file("new.rds", "cache.rds")
    Condition
      Error in `.replace_cache_file()`:
      ! Cannot replace the cache file 'cache.rds'.
