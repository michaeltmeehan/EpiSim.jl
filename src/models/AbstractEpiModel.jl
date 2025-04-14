abstract type AbstractEpiModel end


"""
    update_event_rates!(event_rates, model, ...)

Update the event rates for the given model.
Must be implemented for each subtype of `AbstractEpiModel`.
"""